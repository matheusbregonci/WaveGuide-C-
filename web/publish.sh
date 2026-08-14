#!/usr/bin/env bash
# Assemble web/site/ with ONLY what the browser needs, ready to upload.
#
# The web/ folder mixes three things: what runs (html/css/js/wasm), what builds
# it (bindings.cpp, build.sh) and what tests it (*.mjs). Uploading all of it
# works but publishes build scripts and test harnesses as part of the site, and
# makes it harder to see what the deployment actually consists of.
#
#   ./web/publish.sh      -> docs/
#
# Output is docs/ because GitHub Pages only serves from the repo root or
# from /docs. Putting the site there means Pages works with no CI at all.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$(dirname "$here")/docs"

if [ ! -f "$here/dist/waveguide.wasm" ]; then
    echo "web/dist/waveguide.wasm nao existe — rode o build do WASM antes." >&2
    exit 1
fi

rm -rf "$out"
mkdir -p "$out/dist"
cp "$here"/index.html "$here"/style.css "$out"/
cp "$here"/app.js "$here"/renderer3d.js "$here"/section2d.js "$here"/theory.js "$out"/
cp "$here"/dist/waveguide.js "$here"/dist/waveguide.wasm "$out/dist/"

# Netlify and Cloudflare Pages both read _headers. Two things matter:
#   - .wasm must arrive as application/wasm, or instantiateStreaming refuses it
#     (some hosts still default to application/octet-stream).
#   - the wasm is content that changes when the physics changes, and students
#     would otherwise keep a cached copy of the OLD solver. Short cache on the
#     entry points, long cache only where a stale copy is harmless.
cat > "$out/_headers" <<'HDR'
/dist/*.wasm
  Content-Type: application/wasm
  Cache-Control: public, max-age=300, must-revalidate

/dist/*.js
  Cache-Control: public, max-age=300, must-revalidate

/*.js
  Cache-Control: public, max-age=300, must-revalidate

/index.html
  Cache-Control: public, max-age=60, must-revalidate
HDR

# GitHub Pages runs Jekyll by default, which SKIPS files and folders starting
# with an underscore -- it would silently drop _headers and any future _assets.
touch "$out/.nojekyll"

echo "pronto: $out"
du -sh "$out" 2>/dev/null || true
find "$out" -type f | sed "s|$out/|  |"

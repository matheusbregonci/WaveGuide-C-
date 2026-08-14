# WaveGuide

Visualizador de campos em guias de onda e cavidades — retangulares (TE/TM_mn)
e cilíndricas (TE/TM_nm) — com aplicativo desktop em C++/OpenGL e uma versão
web em WebAssembly.

## O que ele faz

- Campo **E** ou **H** de qualquer modo TE/TM, em guia aberta ou cavidade
- Nuvem de intensidade 3D, cortes XY/ZX/ZY com setas ou linhas de campo
- **Amplitude normalizada por potência transportada**, então os eixos estão em
  V/m e A/m reais — a mesma convenção de uma *wave port* do HFSS (1 W por
  padrão), que é o que torna os dois comparáveis
- Exportação de um projeto LaTeX completo (figuras + tabelas de configuração)
  para uso em artigo ou material didático

A amplitude só é física quando há potência a transportar. Cavidade armazena
energia e modo evanescente não propaga: nesses casos a escala é arbitrária e a
interface diz isso, em vez de imprimir um número em A/m que não significa nada.

## Versão web

`docs/` é um site estático (HTML + WebAssembly, ~165 KB) servido pelo GitHub
Pages. Só o domínio guia/cavidade foi portado — o que deixa de fora Eigen e
OpenMP, e por isso o módulo **não precisa** de `SharedArrayBuffer` nem dos
cabeçalhos COOP/COEP que normalmente acompanham WASM.

A física do navegador é o **mesmo C++** do desktop compilado para WASM, não uma
reimplementação: `web/check.mjs` confere que os números batem.

## Construir

Desktop:

```bash
cmake -S . -B build && cmake --build build -j
```

Para distribuir, desligue o ajuste ao processador local — senão o binário só
roda em CPUs iguais à que compilou:

```bash
cmake -S . -B build-dist -DWGTE_NATIVE_ARCH=OFF && cmake --build build-dist -j
cd build-dist && cpack
```

Web (requer [emsdk](https://emscripten.org/)):

```bash
./web/build.sh      # -> web/dist/waveguide.{js,wasm}
./web/publish.sh    # -> docs/
```

## Testes

Não há framework; são programas que rodam e comparam contra valores conhecidos.

| | o que verifica |
|---|---|
| `web/check.mjs` | os números do WASM contra o binário desktop |
| `web/integration.mjs` | a fronteira JS↔WASM: formatos, strides, normalizações |
| `web/theorycheck.mjs` | as fórmulas exibidas e o realce termo↔região |

## Estrutura

```
include/, src/     modelos analiticos, solvers numericos, render, exportador
  TEmnModel        guia retangular TE/TM_mn
  CylindricalModel guia circular TE/TM_nm
  FieldViz         geometria de visualizacao (sem GL, sem UI) — compartilhada
                   entre desktop, exportador e WASM
  Bessel.hpp       J_n proprio: as funcoes especiais do C++17 nao existem na
                   libc++, que e o que o Emscripten usa
shaders/           GLSL do desktop
web/               ligacoes WASM, front-end e testes
docs/              site publicado (gerado por web/publish.sh)
```

## Referências

As equações seguem Pozar, *Microwave Engineering* (4ª ed.), seção 3.3 para a
guia retangular e 3.4 para a circular. O material de referência usado no
desenvolvimento não é redistribuído aqui.

cmake_minimum_required(VERSION 3.20)

project(WaveguideTEmn
    VERSION 1.0.0
    DESCRIPTION "TEmn Waveguide Field Visualizer (OpenGL + Raytraced Point Cloud)"
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compiler warnings
if(MSVC)
    add_compile_options(/W4 /utf-8)
else()
    add_compile_options(-Wall -Wextra)
endif()

# ------------------------------------------------------------------
# Dependencies (fetched automatically by CMake on first configure)
# ------------------------------------------------------------------
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# --- GLFW ---
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(glfw)

# --- GLM (header only) ---
FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG        1.0.1
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(glm)

# OpenGL: use the system loader (opengl32.dll on Windows).
# We provide our own minimal GL 3.3 function loader in src/gl_loader.cpp so
# we do not need glad (which requires Python + jinja2 at build time).
find_package(OpenGL REQUIRED)

# --- Dear ImGui (UI) ---
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.90.8-docking   # docking branch: enables DockSpace + dockable windows
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(imgui)
if(NOT imgui_POPULATED)
    FetchContent_Populate(imgui)
endif()

add_library(imgui_lib STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui_lib PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui_lib PUBLIC glfw OpenGL::GL)

# --- Eigen (header-only sparse complex linear algebra for the FDFD solver) ---
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)   # header-only: just expose the include dir
endif()

# ------------------------------------------------------------------
# Executable
# ------------------------------------------------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE include ${eigen_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        glfw
        glm::glm
        OpenGL::GL
        imgui_lib
)

# Performance: native tuning + OpenMP so the FDTD step loops run multi-threaded
# (the microstrip/FDTD stepH/stepE are parallel over the outer plane).
#
# -march=native emits instructions for the CPU that COMPILED the binary. That is
# what you want on your own machine and exactly what you must not ship: an older
# processor without AVX-512 (or whatever the build host had) dies with an
# illegal-instruction fault at some random later point. Keep it on for local
# builds, turn it off for anything that leaves this machine.
option(WGTE_NATIVE_ARCH "Tune for the building CPU (-march=native). MUST be OFF for redistributable builds." ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE -O3)
    if(WGTE_NATIVE_ARCH)
        target_compile_options(${PROJECT_NAME} PRIVATE -march=native)
        message(WARNING
            "WGTE_NATIVE_ARCH=ON: this build is tuned to THIS CPU and is not "
            "safe to distribute. Configure with -DWGTE_NATIVE_ARCH=OFF to package.")
    else()
        # x86-64-v2 = SSE4.2/POPCNT: everything since ~2009. Portable in
        # practice while still well ahead of the plain x86-64 baseline.
        target_compile_options(${PROJECT_NAME} PRIVATE -mtune=generic)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|AMD64)")
            target_compile_options(${PROJECT_NAME} PRIVATE -march=x86-64-v2)
        endif()
        message(STATUS "Portable build: generic tuning, safe to redistribute")
    endif()
endif()

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    if(MINGW)
        # -static below covers libstdc++/libgcc, but CMake's OpenMP::OpenMP_CXX
        # puts -fopenmp on the LINK line too, and the driver's implicit -lgomp
        # resolved through the import library -- leaving the exe needing
        # libgomp-1.dll (and libwinpthread-1.dll behind it), the one dependency
        # that stopped the build from being self-contained.
        #
        # Link the archive by path instead. A -Wl,-Bstatic/-Bdynamic pair does
        # NOT work here: the trailing -Bdynamic switches the linker back for
        # everything resolved afterwards, which includes libstdc++, so it
        # quietly undid -static and made the result worse.
        target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
        # Ask the COMPILER for its own libgomp.a rather than searching the
        # system: MSYS2 ships parallel toolchains (mingw64, ucrt64, clang64)
        # whose archives are not interchangeable, and a plain find_library
        # happily returned the ucrt64 copy for a mingw64 build.
        execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgomp.a
            OUTPUT_VARIABLE WGTE_GOMP_STATIC
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        # gcc echoes the bare name back when it cannot find the file.
        if(WGTE_GOMP_STATIC STREQUAL "libgomp.a" OR NOT EXISTS "${WGTE_GOMP_STATIC}")
            set(WGTE_GOMP_STATIC "")
        endif()
        if(WGTE_GOMP_STATIC)
            target_link_libraries(${PROJECT_NAME} PRIVATE ${WGTE_GOMP_STATIC})
            message(STATUS "OpenMP: static ${WGTE_GOMP_STATIC}")
        else()
            target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
            message(WARNING "libgomp.a not found - the binary will need libgomp-1.dll")
        endif()
    else()
        target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    message(STATUS "OpenMP enabled for ${PROJECT_NAME}")
else()
    message(STATUS "OpenMP not found - FDTD step loops will run single-threaded")
endif()

# Statically link the MinGW/GCC C++ runtime (libstdc++, libgcc, libwinpthread)
# so the executable does not depend on the MSYS2/MinGW DLLs at runtime. Without
# this the program runs from a MinGW shell (which has the DLLs on PATH) but
# crashes with an access violation inside libstdc++-6.dll when double-clicked
# from Explorer, because Windows loads an incompatible runtime DLL from PATH.
if(MINGW OR (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND WIN32))
    target_link_options(${PROJECT_NAME} PRIVATE
        -static -static-libgcc -static-libstdc++)
endif()

# Copy shaders next to the executable so the program can load them at runtime.
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/shaders
        $<TARGET_FILE_DIR:${PROJECT_NAME}>/shaders
    COMMENT "Copying shaders to build directory"
)

# ------------------------------------------------------------------
# Install / package
# ------------------------------------------------------------------
# The shaders are read from disk at run time, so they are part of the product,
# not a build artefact: they must land next to the executable in the install
# tree exactly as they do in build/bin.
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/shaders/ DESTINATION shaders)

set(CPACK_PACKAGE_NAME              "WaveguideTEmn")
set(CPACK_PACKAGE_VENDOR            "WaveguideTEmn")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_VERSION           "${PROJECT_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "WaveguideTEmn")
set(CPACK_PACKAGE_EXECUTABLES       "WaveguideTEmn;Waveguide TEmn")
set(CPACK_CREATE_DESKTOP_LINKS      "WaveguideTEmn")
set(CPACK_PACKAGE_FILE_NAME
    "WaveguideTEmn-${PROJECT_VERSION}-win64")

# ZIP always works and needs no extra tooling, so it is the fallback that keeps
# `cpack` useful on a bare machine. NSIS is added only when makensis is actually
# present -- listing a generator whose tool is missing makes cpack fail outright
# instead of producing the packages it CAN build.
set(CPACK_GENERATOR "ZIP")
find_program(MAKENSIS_EXECUTABLE makensis)
if(MAKENSIS_EXECUTABLE)
    list(APPEND CPACK_GENERATOR "NSIS")
    set(CPACK_NSIS_PACKAGE_NAME     "Waveguide TEmn")
    set(CPACK_NSIS_DISPLAY_NAME     "Waveguide TEmn ${PROJECT_VERSION}")
    set(CPACK_NSIS_INSTALLED_ICON_NAME "WaveguideTEmn.exe")
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
    set(CPACK_NSIS_MODIFY_PATH      OFF)
    message(STATUS "NSIS found (${MAKENSIS_EXECUTABLE}) - `cpack` will build an installer")
else()
    message(STATUS "NSIS not found - `cpack` will build the ZIP only. "
                   "Install NSIS to get a .exe installer.")
endif()

if(WGTE_NATIVE_ARCH)
    message(STATUS "NOTE: packaging a WGTE_NATIVE_ARCH=ON build ships a "
                   "CPU-specific binary. Reconfigure with -DWGTE_NATIVE_ARCH=OFF first.")
endif()

include(CPack)

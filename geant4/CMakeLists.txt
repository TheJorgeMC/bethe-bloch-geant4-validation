# ============================================================================
# CMakeLists.txt — Validación de Bethe-Bloch con protones en slab absorbedor
# Requiere Geant4 11.x. El soporte multi-hilo se hereda de la instalación de
# Geant4: si Geant4 fue compilado con GEANT4_BUILD_MULTITHREADED=ON, el
# G4RunManagerFactory de slab.cc creará un run manager MT automáticamente.
# ============================================================================
cmake_minimum_required(VERSION 3.16...3.27)
project(slab)

# ----------------------------------------------------------------------------
# Localizar Geant4 con soporte de UI y visualización (ui_all vis_all).
# ----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Compilar con soporte de UI y visualización" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# ----------------------------------------------------------------------------
# Fuentes y cabeceras
# ----------------------------------------------------------------------------
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(slab slab.cc ${sources} ${headers})
target_link_libraries(slab ${Geant4_LIBRARIES})

# ----------------------------------------------------------------------------
# Copiar las macros al directorio de build para poder ejecutar
#   ./slab macros/run.mac
# directamente desde el directorio de compilación.
# ----------------------------------------------------------------------------
set(SLAB_MACROS
  macros/run.mac
  macros/vis.mac
  macros/energy_scan.mac
  macros/run_point.mac
  macros/generate_energy_scan.py
)
foreach(_macro ${SLAB_MACROS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_macro}
    ${PROJECT_BINARY_DIR}/${_macro}
    COPYONLY
  )
endforeach()

install(TARGETS slab DESTINATION bin)

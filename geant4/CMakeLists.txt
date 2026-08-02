# ============================================================================
# CMakeLists.txt — Bethe-Bloch validation with protons on an absorber slab
# Requires Geant4 11.x. Multi-threading support is inherited from the Geant4
# installation: if Geant4 was built with GEANT4_BUILD_MULTITHREADED=ON, the
# G4RunManagerFactory in slab.cc will create an MT run manager automatically.
# ============================================================================
cmake_minimum_required(VERSION 3.16...3.27)
project(slab)

# ----------------------------------------------------------------------------
# Locate Geant4 with UI and visualization support (ui_all vis_all).
# ----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build with UI and visualization support" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# ----------------------------------------------------------------------------
# Sources and headers
# ----------------------------------------------------------------------------
file(GLOB sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(slab slab.cc ${sources} ${headers})
target_link_libraries(slab ${Geant4_LIBRARIES})

# ----------------------------------------------------------------------------
# Copy the macros into the build directory so that
#   ./slab macros/run.mac
# works directly from the build tree.
# ----------------------------------------------------------------------------
file(GLOB SCAN_MACROS
        RELATIVE ${PROJECT_SOURCE_DIR}
        CONFIGURE_DEPENDS
        "${PROJECT_SOURCE_DIR}/macros/run_scan/run_*MeV.mac"
)

#   Fixed macros plus found by the GLOB scan
set(SLAB_MACROS
        macros/vis.mac
        macros/energy_scan.mac
        macros/run_point.mac
        ${SCAN_MACROS}
)

#run_scan macros loop
foreach(_macro ${SLAB_MACROS})
  configure_file(
          ${PROJECT_SOURCE_DIR}/${_macro}
          ${PROJECT_BINARY_DIR}/${_macro}
          COPYONLY
  )
endforeach()

file(COPY "${PROJECT_SOURCE_DIR}/macros/material_data"
        DESTINATION "${PROJECT_BINARY_DIR}/macros"
)

install(TARGETS slab DESTINATION bin)

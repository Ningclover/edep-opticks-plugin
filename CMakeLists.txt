cmake_minimum_required(VERSION 3.18)
project(edep-simphony-plugin LANGUAGES CXX)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Remember any C++ standard requested on the command line (-DCMAKE_CXX_STANDARD=NN)
# before find_package(EDepSim) below overwrites CMAKE_CXX_STANDARD with ROOT's.
set(RequestedCxxStandard "${CMAKE_CXX_STANDARD}")

# ── Dependencies ──────────────────────────────────────────────────────────────
# Everything is located via find_package / CMAKE_PREFIX_PATH (as set by e.g. a
# Spack environment).  No machine-specific paths or environment variables are
# used.  The imported targets carry their own include directories and transitive
# dependencies:
#   - simphony::* pull in plog, glm, OptiX and the CUDA runtime.
#   - EDepSim::*  pull in the edep-sim headers.
find_package(ROOT REQUIRED COMPONENTS Core Tree RIO)
find_package(Geant4 REQUIRED)
find_package(EDepSim REQUIRED)
# simphony's exported config references glm::glm but omits its own
# find_dependency(glm), so resolve glm here before importing simphony's targets.
find_package(glm REQUIRED)
find_package(simphony REQUIRED)

# ── C++ standard ────────────────────────────────────────────────────────────
# Use the C++ standard given on the command line (-DCMAKE_CXX_STANDARD=NN).  If
# none is given, inherit the standard that ROOT was compiled with: ROOT requires
# that dependent code be built with the same C++ standard as ROOT itself.
set(RootCxxStandard "")
get_target_property(root_features ROOT::Core INTERFACE_COMPILE_FEATURES)
if(root_features)
  foreach(cxxstd IN ITEMS 26 23 20 17 14 11)
    if("cxx_std_${cxxstd}" IN_LIST root_features)
      set(RootCxxStandard ${cxxstd})
      break()
    endif()
  endforeach()
endif()

# Fall back to the standard advertised by ROOTConfig.cmake.
if(NOT RootCxxStandard AND DEFINED ROOT_CXX_STANDARD)
  set(RootCxxStandard ${ROOT_CXX_STANDARD})
endif()

if(RootCxxStandard)
  message(STATUS "ROOT compiled with C++${RootCxxStandard}")
endif()

if(RequestedCxxStandard)
  set(CMAKE_CXX_STANDARD ${RequestedCxxStandard})
  if(RootCxxStandard AND NOT CMAKE_CXX_STANDARD EQUAL RootCxxStandard)
    message(WARNING "Building with C++${CMAKE_CXX_STANDARD} but ROOT requires C++${RootCxxStandard}")
  endif()
elseif(RootCxxStandard)
  set(CMAKE_CXX_STANDARD ${RootCxxStandard})
else()
  message(FATAL_ERROR
    "Could not determine ROOT's C++ standard; set it with -DCMAKE_CXX_STANDARD=NN")
endif()

message(STATUS "Building edep-simphony-plugin with C++${CMAKE_CXX_STANDARD}")

# ── Plugin library ────────────────────────────────────────────────────────────
add_library(edep-simphony-plugin SHARED
    src/SimphonyRunAction.cc
    src/SimphonyEventAction.cc
    src/SimphonyStepAction.cc
    src/SimphonyPhysicsSwap.cc
    src/SimphonyScintProcess.cc
    src/SimphonyCpuPhotonTracker.cc
    src/plugin_entry.cc
)

target_compile_definitions(edep-simphony-plugin PRIVATE
    GLM_ENABLE_EXPERIMENTAL OPTICKS_U4 STANDALONE)

target_include_directories(edep-simphony-plugin PRIVATE src)

target_link_libraries(edep-simphony-plugin PRIVATE
    # Geant4 + ROOT
    ${Geant4_LIBRARIES}
    ROOT::Core ROOT::Tree ROOT::RIO
    # edep-sim
    EDepSim::edepsim
    EDepSim::edepsim_io
    # simphony (formerly eic-opticks); order matches static-dep order.
    simphony::G4CX
    simphony::U4
    simphony::CSGOptiX
    simphony::QUDARap
    simphony::CSG
    simphony::SysRap
)

install(TARGETS edep-simphony-plugin LIBRARY DESTINATION lib)

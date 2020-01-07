CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)


# ====================================
# Collect Source Files
# ====================================
FILE(GLOB PhaseTransition_SOURCE *.cpp)

# ====================================
# Build the libraries
# ====================================
# libphasetransition.a
ADD_LIBRARY(PhaseTransition_static STATIC ${PhaseTransition_SOURCE})
SET_TARGET_PROPERTIES(PhaseTransition_static PROPERTIES LINKER_LANGUAGE CXX
    VERSION ${PhaseTransition_VERSION} OUTPUT_NAME phasetransition
)
TARGET_LINK_LIBRARIES(PhaseTransition_static ${GSL_LIBRARIES})

# libphasetransition.so
ADD_LIBRARY(PhaseTransition_shared SHARED ${PhaseTransition_SOURCE})
SET_TARGET_PROPERTIES(PhaseTransition_shared PROPERTIES LINKER_LANGUAGE CXX
    VERSION ${PhaseTransition_VERSION} OUTPUT_NAME phasetransition
)
TARGET_LINK_LIBRARIES(PhaseTransition_shared ${GSL_LIBRARIES})

SET(PhaseTransition_LIBRARIES PhaseTransition_shared)

# SET(PhaseTransition_INSTALL_LIBS PhaseTransition_static PhaseTransition_shared)
# INSTALL(TARGETS ${PhaseTransition_INSTALL_LIBS} DESTINATION )

add_subdirectory(unit_test)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)


add_executable(HermiteInterpolationTest check_HI.cpp)
add_executable(PhaseTest check_Phase.cpp)
add_executable(TraceTest check_traceMin.cpp)
add_executable(TraceMultiTest check_traceMultiMin.cpp)
add_executable(TcTest check_findTc.cpp)
add_executable(MemTest check_MemeryUsage.cpp)
add_executable(LocalTest check_ApproxLocalMin.cpp)
add_executable(Tunnel1D check_Tunneling1D.cpp)
add_executable(functiontest check_functional.cpp)


target_link_libraries(HermiteInterpolationTest ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})
target_link_libraries(PhaseTest ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})
target_link_libraries(TraceTest ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})
target_link_libraries(TraceMultiTest ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})
target_link_libraries(TcTest ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})
target_link_libraries(MemTest ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})
target_link_libraries(LocalTest ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})
target_link_libraries(Tunnel1D ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${PhaseTransition_LIBRARIES})

# install(TARGETS HermiteInterpolationTest DESTINATION ${PROJECT_BINARY_DIR}/bin)
# install(TARGETS PhaseTest DESTINATION ${PROJECT_BINARY_DIR}/bin)
# install(TARGETS TraceTest DESTINATION ${PROJECT_BINARY_DIR}/bin)
# install(TARGETS TraceMultiTest DESTINATION ${PROJECT_BINARY_DIR}/bin)
# install(TARGETS TcTest DESTINATION ${PROJECT_BINARY_DIR}/bin)
# install(TARGETS LocalTest DESTINATION ${PROJECT_BINARY_DIR}/bin)
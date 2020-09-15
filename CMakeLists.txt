cmake_minimum_required (VERSION 2.8)

project(MultilevelExahypeProject)

find_package(MUQ REQUIRED)

include_directories(${MUQ_INCLUDE_DIRS})

message("MUQ_INCLUDE_DIRS: ${MUQ_INCLUDE_DIRS}")
message("MUQ_LIBRARIES: ${MUQ_LIBRARIES}")

set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})
set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})


file(READ "${EXAHYPE_PATH}/exahype_link" EXAHYPE_LINKS)
string(STRIP ${EXAHYPE_LINKS} EXAHYPE_LINKS)

configure_file(exahype.log-filter exahype.log-filter COPYONLY)

add_custom_target(link_target ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/Output ${CMAKE_CURRENT_BINARY_DIR}/Output)

add_executable(SLMCMC SLMCMC.cpp)
target_link_libraries(SLMCMC ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
target_include_directories(SLMCMC PRIVATE ${EXAHYPE_PATH})
set_target_properties(SLMCMC PROPERTIES LINK_FLAGS ${EXAHYPE_LINKS})

add_executable(MLMCMC MLMCMC.cpp)
target_link_libraries(MLMCMC ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
target_include_directories(MLMCMC PRIVATE ${EXAHYPE_PATH})
set_target_properties(MLMCMC PROPERTIES LINK_FLAGS ${EXAHYPE_LINKS})

add_executable(FullParallelMLMCMC FullParallelMLMCMC.cpp)
target_link_libraries(FullParallelMLMCMC ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
target_include_directories(FullParallelMLMCMC PRIVATE ${EXAHYPE_PATH})
set_target_properties(FullParallelMLMCMC PROPERTIES LINK_FLAGS ${EXAHYPE_LINKS})

add_executable(MCMC_LOH1 MCMC_LOH1.cpp)
target_link_libraries(MCMC_LOH1 ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
target_include_directories(MCMC_LOH1 PRIVATE ${EXAHYPE_PATH})
set_target_properties(MCMC_LOH1 PROPERTIES LINK_FLAGS ${EXAHYPE_LINKS})

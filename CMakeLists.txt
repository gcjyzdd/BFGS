#
# Author: Changjie Guan (changjieguan@gmail.com)
#
#
#BSD 2-Clause License
#
#Copyright (c) 2018, Changjie
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#* Redistributions of source code must retain the above copyright notice, this
#  list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
project(bfgs_mpc)

cmake_minimum_required (VERSION 3.5)


#Flags for compiler
IF(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
    set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -mavx2 -mfma -march=native -pg")
ENDIF()


# Set output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(LIB_NAME bfgs_lib)

include_directories(/usr/local/include
					${PROJECT_SOURCE_DIR}/src/Eigen-3.3
					${PROJECT_SOURCE_DIR}/src/bfgs
					${PROJECT_SOURCE_DIR}/src/mpc
                    ${PROJECT_SOURCE_DIR}/src/phr
)

link_directories(/usr/local/lib)

# OpenMP
OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

add_subdirectory(src/bfgs)
add_subdirectory(src/mpc)
add_subdirectory(src/phr)

set(SRC src/main.cpp)
set(SRC_Test_BFGS src/test_struct_function.cpp)
set(SRC_Test_PHR src/test_phr.cpp)
set(SRC_Test_PHR_MPC src/test-mpc-phr/test_mpc_phr.cpp)

add_executable(bfgs_test ${SRC_Test_BFGS})
add_executable(mpc_test ${SRC})
add_executable(phr_test ${SRC_Test_PHR})
add_executable(phr_mpc_test ${SRC_Test_PHR_MPC})

target_link_libraries(bfgs_test bfgs python2.7)
target_link_libraries(mpc_test bfgs mpc python2.7)
target_link_libraries(phr_test bfgs phr)
target_link_libraries(phr_mpc_test bfgs phr)

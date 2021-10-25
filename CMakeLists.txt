cmake_minimum_required(VERSION 2.8)

project (HeartNext)
set(CMAKE_BUILD_TYPE Debug)
set(linux_flags -g  -std=c++0x  -fpic -Wall )
add_definitions(${linux_flags})
FIND_PACKAGE(PkgConfig)
pkg_check_modules(lua REQUIRED lua5.1)
include_directories(${lua_INCLUDE_DIRS})

FIND_PACKAGE(CUDA)
set(CUDA_NVCC_FLAGS  ${CUDA_NVCC_FLAGS} -maxrregcount=50 -arch=compute_20 -Xcompiler -fpic  -Xcompiler -Wall )
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# the vector
cuda_compile(baseCUDA src/base.cu)
add_library(vector SHARED src/Lua.cc src/VectorWrapper.cpp ${baseCUDA})
target_link_libraries(vector lua5.1 ${CUDA_LIBRARIES})

# the 06model
cuda_compile(06model_cu src/06Model.cu  )
add_library(06model SHARED src/06Model.cpp src/Lua.cc ${06model_cu})
target_link_libraries(06model lua5.1 ${CUDA_LIBRARIES})
# the block kr cal model
cuda_compile(block_kr_cal_cu src/BlockKr_Cal.cu)
add_library(block_kr_cal SHARED src/BlockKr_Cal.cpp src/Lua.cc ${block_kr_cal_cu})
# the purkinje model
cuda_compile(purkinjeModel_cu src/PurkinjeModel.cu)
add_library(purkinjeModel SHARED src/PurkinjeModel.cpp ${purkinjeModel_cu})
target_link_libraries(purkinjeModel lua5.1 ${CUDA_LIBRARIES})

# the test programe
add_executable(test src/Lua.cc src/test.cpp)
target_link_libraries(test lua5.1 vector 06model purkinjeModel block_kr_cal)
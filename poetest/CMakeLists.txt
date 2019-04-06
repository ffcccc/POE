cmake_minimum_required(VERSION 2.8.9)
project (poetest)
 
#For the shared library:
#set ( PROJECT_LINK_LIBS poelib )
#link_directories( ~/exploringBB/extras/cmake/studentlib_shared/build )
 
#For the static library:
set ( PROJECT_LINK_LIBS poelib )
link_directories( ../build )

set ( PROJECT_INC_BOOST C:/Users/735011/workspace/3rdparty/boost )
set ( PROJECT_INC_EIGEN C:/Users/735011/workspace/3rdparty/eigen3 )
set ( PROJECT_INC_UTILS C:/Users/735011/workspace/KCL/EigenUtils )
set ( PROJECT_INC_CSV C:/Users/735011/workspace/3rdparty/fast-csv-parser-master )
set ( PROJECT_INCS C:/Users/735011/workspace/3rdparty )

include_directories(${PROJECT_INCS})
include_directories(${PROJECT_INC_BOOST} ${PROJECT_INC_EIGEN} ${PROJECT_INC_UTILS} ${PROJECT_INC_CSV})
include_directories(../poelib/include)

add_executable(poetest test.cpp)
target_link_libraries(poetest ${PROJECT_LINK_LIBS} )
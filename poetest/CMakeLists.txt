cmake_minimum_required(VERSION 2.8.9)
project (poetest)
 
#For the shared library:
#set ( PROJECT_LINK_LIBS poelib.so )
#link_directories( ~/exploringBB/extras/cmake/studentlib_shared/build )
 
#For the static library:
set ( PROJECT_LINK_LIBS poelib.a )
link_directories( ../poelib/build )
 
include_directories(../poelib/include)
 
add_executable(poetest test.cpp)
target_link_libraries(poetest ${PROJECT_LINK_LIBS} )
cmake_minimum_required(VERSION 2.8.9)
project(poelib)
set(CMAKE_BUILD_TYPE Release)
 
#Bring the headers, such as Student.h into the project
set ( PROJECT_INC_BOOST C:/Users/735011/workspace/3rdparty/boost )
set ( PROJECT_INC_EIGEN C:/Users/735011/workspace/3rdparty/eigen3 )
set ( PROJECT_INC_UTILS C:/Users/735011/workspace/KCL/EigenUtils )
set ( PROJECT_INC_CSV C:/Users/735011/workspace/3rdparty/fast-csv-parser-master )
set ( PROJECT_INCS C:/Users/735011/workspace/3rdparty )

include_directories(${PROJECT_INCS})
include_directories(${PROJECT_INC_BOOST} ${PROJECT_INC_EIGEN} ${PROJECT_INC_UTILS} ${PROJECT_INC_CSV})
include_directories(include)
 
#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")
file(GLOB NETSTAT "src/dcdflib/*.c")
#Generate the static library from the sources
add_library(poelib STATIC ${SOURCES} ${NETSTAT})
 
#Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply
install(TARGETS poelib DESTINATION lib)
cmake_minimum_required(VERSION 3.8)
project(p VERSION 1.0 LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 11)

set(GAMB_LIB /cs/labs/dina/dina/libs/gamb)
set(DOCK_LIB /cs/labs/dina/dina/libs/DockingLib)
set(EM_LIB /cs/labs/dina/dina/libs/EMlib)
set(LEDA_LIB /cs/labs/dina/dina/libs/leda/)
set(LEDA_INCL /cs/labs/dina/dina/libs/leda/incl/)



set(SOURCE_FILES
        BestK.h
        TransformationSearch.cpp
        TransformationSearch.h
        FFTCommon.cpp
        FFTCommon.h
        main.cpp)

include_directories(${LEDA_INCL} ${GAMB_LIB} ${DOCK_LIB} ${EM_LIB})
link_directories(${GAMB_LIB} ${DOCK_LIB} ${LEDA_LIB} ${EM_LIB})


add_executable(p ${SOURCE_FILES})

target_link_libraries(p emLib dockingLib gamb++ stdc++ leda X11 m -lfftw3_threads -lfftw3 -lm -pthread)
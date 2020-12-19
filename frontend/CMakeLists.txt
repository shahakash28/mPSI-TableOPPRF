
#project(fronend)




#############################################
#               Build fronend.exe           #
#############################################

file(GLOB_RECURSE SRC_FRONTEND ${CMAKE_SOURCE_DIR}/frontend/*.cpp)
include_directories(${CMAKE_SOURCE_DIR}/frontend/)

add_executable(frontend.exe ${SRC_FRONTEND})

target_include_directories(frontend.exe PUBLIC $ENV{HOME}/libscapi/build/NTL/include $ENV{HOME} $ENV{HOME}/libscapi/install/include
        $ENV{HOME}/libscapi/install/include/OTExtensionBristol)

#############################################
#            Link our libraries             #
#############################################

#target_link_libraries(frontend.exe libPSI_Tests)
target_link_libraries(frontend.exe libOPRF)
target_link_libraries(frontend.exe libOTe)
target_link_libraries(frontend.exe gmp)
target_link_libraries(frontend.exe $ENV{HOME}/libscapi/libscapi.a)
target_link_libraries(frontend.exe $ENV{HOME}/libscapi/build/NTL/src/ntl.a)
#        $ENV{HOME}/libscapi/libscapi.a gmp OTExtensionBristol libsimpleot.a
target_link_libraries(frontend.exe gmp)
target_link_libraries(frontend.exe boost_system boost_thread pthread crypto dl ssl z)

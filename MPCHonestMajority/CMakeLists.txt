cmake_minimum_required(VERSION 2.8)
project(MPCHonestMajority)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -maes -mpclmul -mbmi2 -msse4.1 -Wno-ignored-attributes")


INCLUDE_DIRECTORIES($ENV{HOME}/libscapi/build/NTL/include $ENV{HOME} $ENV{HOME}/libscapi/install/include
        $ENV{HOME}/libscapi/install/include/OTExtensionBristol)

set(SOURCE_FILES ProtocolParty.cpp ProtocolParty2.cpp TemplateField.cpp ZpKaratsubaElement.cpp main.cpp)


link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib )


add_executable(MPCHonestMajority ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(MPCHonestMajority $ENV{HOME}/libscapi/build/NTL/src/ntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp OTExtensionBristol libsimpleot.a
        boost_system boost_thread pthread crypto dl ssl z)

include_guard()

file(GLOB SOURCES src/*.cpp)
file(GLOB INCLUDES includes/*.h)

add_executable(TestSuite ${SOURCES})

find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
target_include_directories (TestSuite 
    PUBLIC 
        ${LINEAR_SOURCE_DIR}/src
        ${LINEAR_SOURCE_DIR}/includes
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    )

add_definitions (-DBOOST_TEST_DYN_LINK)

target_compile_options(TestSuite PRIVATE -W -Wall -fno-elide-constructors)

target_link_libraries (TestSuite
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_test (TestSuite TestSuite)
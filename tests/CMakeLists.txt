find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories (${LINEAR_SOURCE_DIR}/src
                     ${LINEAR_SOURCE_DIR}/includes
                     ${Boost_FILESYSTEM_LIBRARY}
                     ${Boost_SYSTEM_LIBRARY}
                     ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                    )

add_definitions (-DBOOST_TEST_DYN_LINK)

file(GLOB SOURCES src/*.cpp)

add_executable(TestSuite ${SOURCES})

target_link_libraries (TestSuite
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )
target_compile_options (TestSuite
                        PUBLIC -Wall)


add_test (TestSuite TestSuite)
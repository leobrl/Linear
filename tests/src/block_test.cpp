#define BOOST_TEST_MODULE example
#include <boost/test/included/unit_test.hpp>
#include <block.hpp>

BOOST_AUTO_TEST_SUITE( test_suite_1 )

BOOST_AUTO_TEST_CASE( test_block_creation_1,
                      * boost::unit_test::description("description test 1") )
{
  BOOST_TEST( true /* test assertion */ );
}

BOOST_AUTO_TEST_SUITE_END()
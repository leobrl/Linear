#define BOOST_TEST_MODULE example
#include <boost/test/included/unit_test.hpp>
#include <Linear.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE( test_suite_1 )

BOOST_AUTO_TEST_CASE( test_block_creation,
                      * boost::unit_test::description("Tests instantiation of memory block.") )
{
  const size_t n {10};
  auto block{ linear::Block<int>(n) };

  BOOST_TEST( true );
}

BOOST_AUTO_TEST_CASE( test_block_fill,
                      * boost::unit_test::description("Tests initialization of memory block\
                      * to const values.") )
{
  const size_t n {10};
  auto value {5.2}; 
  auto block{ linear::Block<decltype(value)>(n) };
  block.fill(value);

  auto pass {true};
  for (size_t idx=0; idx<n; ++idx){
    pass = block(idx) == value ? pass : pass & false;
  }

  BOOST_TEST( pass );
}

BOOST_AUTO_TEST_CASE( test_block_assignment,
                      * boost::unit_test::description("Tests memory block assignment.") )
{
  const size_t n {10};
  auto one {1}; 
  auto two {2};

  auto rhs{ linear::Block<decltype(one)>(n) };
  rhs.fill(one);

  linear::Block<decltype(one)> lhs {};

  lhs = rhs;
  rhs.fill(two);

  auto pass {true};
  for (size_t idx=0; idx<n; ++idx){
    pass = lhs(idx) == one ? pass : pass & false;
    pass = rhs(idx) == two ? pass : pass & false;
  }

  BOOST_TEST( pass );
}

BOOST_AUTO_TEST_SUITE_END()
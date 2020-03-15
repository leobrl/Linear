#define BOOST_TEST_MODULE example
#include <boost/test/included/unit_test.hpp>
#include <Linear.hpp>
#include <iostream>

typedef linear::Block<int> IntBlock; 

struct IntBlockFixture {
  IntBlockFixture() : n(10) { int_block = IntBlock(n); }
  ~IntBlockFixture() = default;

  bool assert_is_constant(const IntBlock& block, int value){
    auto pass {true};
    for (size_t idx=0; idx<n; ++idx){
      pass = block(idx) == value ? pass : pass & false;
    }

    return(pass);
  }

  size_t n;
  IntBlock int_block;
};

BOOST_FIXTURE_TEST_SUITE(test_block, IntBlockFixture)

BOOST_AUTO_TEST_CASE( test_block_fill,
                      * boost::unit_test::description("Tests write/read of memory block\
                      * to const values.") )
{
  int value {5};
  int_block.fill(value);

  BOOST_TEST( assert_is_constant(int_block, value) );
}

BOOST_AUTO_TEST_CASE( test_block_assignment,
                      * boost::unit_test::description("Tests block assignment operator.") )
{
  auto one {1}; 
  auto two {2};

  int_block.fill(one);

  IntBlock lhs {};

  lhs = int_block;
  int_block.fill(two);

  BOOST_TEST( assert_is_constant(lhs, one) & assert_is_constant(int_block, two) );
}

BOOST_AUTO_TEST_CASE( test_block_copy_ctor,
                      * boost::unit_test::description("Tests block copy constructor.") )
{
  auto one {1}; 
  auto two {2};

  int_block.fill(one);

  auto lhs {IntBlock(int_block)};
  int_block.fill(two);

  BOOST_TEST( assert_is_constant(lhs, one) & assert_is_constant(int_block, two) );
}

BOOST_AUTO_TEST_CASE( test_move_ctor,
                      * boost::unit_test::description("Tests block move constructor.") )
{
  auto one {1};
  
  auto rhs {int_block};
  rhs.fill(one);

  auto lhs = std::move(rhs);

  BOOST_TEST( assert_is_constant(lhs, one));
}

BOOST_AUTO_TEST_SUITE_END()
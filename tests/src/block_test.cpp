#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <Linear.hpp>
#include <iostream>

typedef boost::mpl::list<int, double, float, char> test_types;
typedef linear::Block<int> IntBlock; 

bool assert_is_equal(const IntBlock& block, std::vector<int> v){
  bool pass = v.size() == block.size();
  
  auto it = v.begin();
  for (size_t idx=0; pass & (idx<block.size()) & (it != v.end()); ++idx, ++it){
    pass = pass & (block[idx] == *it);
  }
  
  return pass;
}

template<typename T>
bool assert_is_equal_to_const(const linear::Block<T>& block, T value){
  auto pass {true};
  for (size_t idx=0; pass & (idx<block.size()) ; ++idx){
    pass = pass & (block[idx] == value);
  }
  
  return(pass);
}

struct IntBlockFixture {
  IntBlockFixture() : n(10) { int_block = IntBlock(n); }
  ~IntBlockFixture() = default;

  size_t n;
  IntBlock int_block;
};

BOOST_FIXTURE_TEST_SUITE(block_basic_tests, IntBlockFixture)

BOOST_AUTO_TEST_CASE(test_vector_size,
                      * boost::unit_test::description("Tests constructor from vector.") )
{
  BOOST_TEST(int_block.size() == 10);
}

BOOST_AUTO_TEST_CASE( test_fixture_assert_is_equal_to_const,
                      * boost::unit_test::description("Tests assert is equal to const in fixture." ))
{
  int zero{0};
  int one{1};
  BOOST_TEST((assert_is_equal_to_const(int_block, zero)) & 
              !assert_is_equal_to_const(int_block, one));
}

BOOST_AUTO_TEST_CASE( test_fixture_assert_is_equal_to_vec,
                      * boost::unit_test::description("Tests assert is equal to const in fixture." ))
{
  std::vector<int> zeros {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<int> ones {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<int> shorter_vec {0, 0, 0};
  std::vector<int> longer_vec {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  BOOST_TEST(assert_is_equal(int_block, zeros) & 
              !assert_is_equal(int_block, ones) & 
              !assert_is_equal(int_block, shorter_vec),
              !assert_is_equal(int_block, longer_vec));
}

BOOST_AUTO_TEST_CASE( test_block_fill,
                      * boost::unit_test::description("Tests write/read of memory block\
                      * to const values.") )
{
  int one {1};
  int_block.fill(one);

  BOOST_TEST( assert_is_equal_to_const(int_block, one));
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

  BOOST_TEST( assert_is_equal_to_const(lhs, one) & 
              assert_is_equal_to_const(int_block, two));
}

BOOST_AUTO_TEST_CASE( test_block_copy_ctor,
                      * boost::unit_test::description("Tests block copy constructor.") )
{
  auto one {1}; 
  auto two {2};

  int_block.fill(one);

  auto lhs {IntBlock(int_block)};
  int_block.fill(two);

  BOOST_TEST( assert_is_equal_to_const(lhs, one) & assert_is_equal_to_const(int_block, two) );
}

BOOST_AUTO_TEST_CASE( test_move_ctor,
                      * boost::unit_test::description("Tests block move constructor.") )
{
  auto one {1};
  
  auto rhs {int_block};
  rhs.fill(one);

  auto lhs = std::move(rhs);

  BOOST_TEST( assert_is_equal_to_const(lhs, one));
}

BOOST_AUTO_TEST_CASE(test_vector_ctor,
                      * boost::unit_test::description("Tests constructor from vector.") )
{
  std::vector<int> v {1,2,3,4,5};
  linear::Block<int> b = linear::Block<int>(v);

  BOOST_TEST(assert_is_equal(b, v));
}

BOOST_AUTO_TEST_CASE(test_ostream,
                      * boost::unit_test::description("Tests ostream.") )
{
  int_block.fill(1);

  std::string expected = "1 1 1 1 1 1 1 1 1 1";

  std::ostringstream stream;
  stream << int_block;

  std::string ostr = stream.str();
  auto check_1 = (ostr.length() == expected.length());
  auto check_2 = (std::equal(ostr.begin(), ostr.end(), expected.begin()));
  
  BOOST_TEST( check_1 & check_2 );
}

BOOST_AUTO_TEST_SUITE_END()

//

BOOST_AUTO_TEST_SUITE(block_template_tests)

BOOST_AUTO_TEST_CASE_TEMPLATE( test_template_block_ctor, T, test_types )
{
  linear::natural n{10};
  auto block {linear::Block<T>(n)};
  T value {};

  BOOST_TEST(assert_is_equal_to_const(block, value));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_template_copy_ctor, T, test_types )
{
  linear::natural n{10};
  auto block {linear::Block<T>(n)};
  auto block_copy {block};

  BOOST_TEST(!(&block[0] == &block_copy[0]));

}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_alignment, T, test_types )
{
  linear::natural n{10};
  auto block {linear::Block<T>(n)};
  auto p{block.pfront()};
  
  BOOST_TEST(((unsigned long)p & 15) == 0);  
}

BOOST_AUTO_TEST_SUITE_END()
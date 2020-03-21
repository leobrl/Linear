#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <Linear.hpp>	
#include <iostream>

typedef boost::mpl::list<int, double, float, char> test_types;

BOOST_AUTO_TEST_SUITE(matrix_basic_test)

BOOST_AUTO_TEST_CASE(test_matrix_ostream,
                      * boost::unit_test::description("Tests matrix ostream operator") )
{
	std::string expected = "0 0 0\n0 0 0\n0 0 0";

	linear::natural n_row {3};
	linear::natural n_col {3};
	auto m {linear::Matrix<int>(n_row, n_col)};

	std::ostringstream stream;
	stream << m;
	std::string str = stream.str();
	
	auto check_1 = (str.length() == expected.length());
  	auto check_2 = (std::equal(str.begin(), str.end(), expected.begin()));

	BOOST_TEST(check_1 & check_2);
}

BOOST_AUTO_TEST_CASE(test_vector_ctor,
                      * boost::unit_test::description("Tests diagonal matrix constructor") )
{
	
	std::string expected = "1 0 0\n0 1 0\n0 0 1";
	std::vector<int> diag {1, 1, 1};

	auto m {linear::Matrix<int>(diag)};

	std::ostringstream stream;
	stream << m;
	std::string str = stream.str();

	auto check_1 = (str.length() == expected.length());
  	auto check_2 = (std::equal(str.begin(), str.end(), expected.begin()));

	BOOST_TEST(check_1 & check_2);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(matrix_template_tests)

BOOST_AUTO_TEST_CASE_TEMPLATE(test_base_ctor, T, test_types){
	
	auto mat = linear::Matrix<T>();
	auto value = linear::natural();

	auto check_1 = mat.n_row == value;
	auto check_2 = mat.n_col == value;

	BOOST_TEST(check_1 & check_2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_copy_constructor, T, test_types){
	
	linear::natural n {5};
	auto mat {linear::Matrix<T>(n, n)};
	auto mat_copy {mat};

	BOOST_TEST(not (&mat(0,0) == &mat_copy(0,0)));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_move_op, T, test_types){
	
	linear::natural n {5};
	auto mat {linear::Matrix<T>(n, n)};
	auto ptr = &(mat(0,0));
	
	auto mat_copy = std::move(mat);

	BOOST_TEST(&(mat_copy(0,0)) == ptr);
}

BOOST_AUTO_TEST_SUITE_END()

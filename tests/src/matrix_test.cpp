#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <Linear.hpp>	
#include <iostream>

using test_types = boost::mpl::list<int, double, float, char>;

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

BOOST_AUTO_TEST_CASE(test_diag_ctor,
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

BOOST_AUTO_TEST_CASE(test_vector_ctor,
                      * boost::unit_test::description("Tests matrix constructor from vector of elements") )
{
	
	std::string expected = "1 2 3 4\n5 6 7 8\n";
	std::vector<int> v {1, 2, 3, 4, 5, 6, 7, 8};

	auto m {linear::Matrix<int>(2, 4, v)};

	std::ostringstream stream;
	stream << m;
	std::string str = stream.str();

	auto check_1 = (str.length() == expected.length());
  	auto check_2 = (std::equal(str.begin(), str.end(), expected.begin()));

	BOOST_TEST(check_1 & check_2);
}

BOOST_AUTO_TEST_CASE(test_matrix_transpose,
					* boost::unit_test::description("Tests matrix transpose."))
{
	std::vector<int> m {1, 2, 3, 4, 5, 6}; 
	auto rhs {linear::Matrix<int>(2, 3, m)};

	// 1 2 3 => 1 4
	// 4 5 6	2 5
	//			3 6

	rhs.transpose();
	
	//std::cout << rhs << std::endl;

	std::vector<int> expected {1, 4, 2, 5, 3, 6}; 

	int i {0};
	auto is_correct {true};
	for (auto it = rhs.begin(); it != rhs.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	

	BOOST_TEST(is_correct);

}

BOOST_AUTO_TEST_CASE(test_inplace_matrix_multiplication,
					* boost::unit_test::description("Tests matrix in place matrix multiplication."))
{
	std::vector<int> m {1,2,3,4}; 
	auto rhs {linear::Matrix<int>(2, 2, m)};
	auto lhs {linear::Matrix<int>(2, 2, m)};

	//	1 2	 *	1 2		=	7 	10
	//	3 4		3 4			15 	22

	rhs *= lhs;
	std::vector<int> expected {7, 10, 15, 22}; 

	int i {0};
	auto is_correct {true};
	for (auto it = rhs.begin(); it != rhs.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	

	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(test_matrix_multiplication,
					* boost::unit_test::description("Tests matrix matrix multiplication."))
{
	std::vector<int> m_rhs {1, 2, 3, 4, 5, 6};
	std::vector<int> m_lhs {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}; 
	auto rhs {linear::Matrix<int>(2, 3, m_rhs)};
	auto lhs {linear::Matrix<int>(3, 4, m_lhs)};

	//	1 2 3	 *	1  2  3  4		=	38	44	 50	 56
	//	4 5	6		5  6  7  8			83 	98	113	128
	//				9 10 11 12
	
	auto res = rhs * lhs;
	std::vector<int> expected {38, 44, 50, 56, 83, 98, 113, 128}; 

	int i {0};
	
	auto is_correct {true};
	for (auto it = res.begin(); it != res.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	
	
	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(test_matrix_inplace_matrix_sum,
					* boost::unit_test::description("Tests matrix in place matrix sum."))
{
	std::vector<int> m_rhs {1, 2, 3, 4, 5, 6};
	std::vector<int> m_lhs {0, 1, 2, 3, 4, 5}; 
	auto rhs {linear::Matrix<int>(2, 3, m_rhs)};
	auto lhs {linear::Matrix<int>(2, 3, m_lhs)};
	
	rhs += lhs;
	std::vector<int> expected {1, 3, 5, 7, 9, 11}; 

	int i {0};
	
	auto is_correct {true};
	for (auto it = rhs.begin(); it != rhs.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	
	
	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(test_matrix_inplace_matrix_subtraction,
					* boost::unit_test::description("Tests matrix in place matrix subtraction."))
{
	std::vector<int> m_rhs {1, 2, 3, 4, 5, 6};
	std::vector<int> m_lhs {0, 1, 2, 3, 4, 5}; 
	auto rhs {linear::Matrix<int>(2, 3, m_rhs)};
	auto lhs {linear::Matrix<int>(2, 3, m_lhs)};
	
	rhs -= lhs;
	std::vector<int> expected {1, 1, 1, 1, 1, 1}; 

	int i {0};
	
	auto is_correct {true};
	for (auto it = rhs.begin(); it != rhs.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	
	
	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(test_matrix_instream_operator,
					* boost::unit_test::description("Tests matrix instream operator."))
{
	std::vector<int> vec {1, 2, 3, 4, 5, 6};
	linear::Matrix<int> m(3, 2);

	m >> vec;

	bool is_correct{true};
	auto vec_it = vec.begin();
	for(auto it = m.begin(); it != m.end(); ++it, ++vec_it){
		is_correct &= *it == *vec_it;
	}

	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(test_matrix_instream_operator_exception,
					* boost::unit_test::description("Tests matrix instream operator exception."))
{
	std::vector<int> vec {1, 2, 3, 4, 5, 6, 7};
	linear::Matrix<int> m(3, 2);

	BOOST_CHECK_THROW(m >> vec, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_invalid_vector_ctor,
					 * boost::unit_test::description("Tests matrix constructor from invalid vector of elements")){

	std::vector<int> v {1, 2, 3, 4, 5, 6, 7, 8};
	BOOST_CHECK_THROW(linear::Matrix<int>(2, 5, v), std::invalid_argument);
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(matrix_template_tests)

BOOST_AUTO_TEST_CASE_TEMPLATE(test_base_ctor, T, test_types){
	
	auto mat = linear::Matrix<T>();
	auto value = linear::natural();

	auto check_1 = mat.nrow() == value;
	auto check_2 = mat.ncol() == value;

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

BOOST_AUTO_TEST_CASE_TEMPLATE(test_row_iterator, T, test_types){
	auto vec = std::vector<T>(10);
	
	for(size_t e = 0; e<vec.size(); ++e){
		vec[e] = static_cast<T>(e);
	}

	auto mat = linear::Matrix<T>(2, 5, vec);

	size_t i = 0;

	auto assert_increment {true};
	for (auto it = mat.begin(); it != mat.end(); ++it, ++i){
		assert_increment &= *it == vec[i];
	}
	
	auto assert_decrement {true};
	for (auto it = (--mat.end()); it != mat.begin(); --it, --i){
		assert_decrement &= *it == vec[i-1];
	}

	BOOST_TEST(assert_increment & assert_decrement);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_row_iterator_post_increment, T, test_types){
	auto vec = std::vector<T>(10);
	
	for(size_t e = 0; e<vec.size(); ++e){
		vec[e] = static_cast<T>(e);
	}

	auto mat = linear::Matrix<T>(2, 5, vec);
	
	size_t i = 0;

	auto assert_increment {true};
	for (auto it = mat.begin(); it != mat.end(); it++, ++i){
		assert_increment &= *it == vec[i];
	}
	
	auto assert_decrement {true};
	for (auto it = (--mat.end()); it != mat.begin(); it--, --i){
		assert_decrement &= *(it) == vec[i-1];
	}

	BOOST_TEST(assert_increment & assert_decrement);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_col_iterator, T, test_types){

	auto vec = std::vector<T>(10);
	for(size_t e = 0; e<vec.size(); ++e){
		vec[e] = static_cast<T>(e);
	}

	linear::natural n_col{5};
	linear::natural n_row{2};
	auto mat = linear::Matrix<T>(n_row, n_col, vec);

	auto expected = std::vector<T>();
	for (linear::natural c = 0; c < n_col; ++c){
		for (linear::natural r = 0; r < n_row; ++r){
			expected.push_back(vec[c + r*n_col]); 
		}	
	}
		
	size_t i = 0;
	auto assert_increment {true};
	for (auto it = mat.col_begin(); it != mat.col_end(); ++it, ++i){
		assert_increment &= *it == expected[i];
	}

	auto assert_decrement {true};
	for (auto it = (--mat.col_end()); it != mat.col_begin(); --it, --i){
		assert_decrement &= *(it) == expected[i-1];
	}

	BOOST_TEST(assert_increment);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_col_iterator_post_increment, T, test_types){

	auto vec = std::vector<T>(10);
	for(size_t e = 0; e<vec.size(); ++e){
		vec[e] = static_cast<T>(e);
	}

	linear::natural n_col{5};
	linear::natural n_row{2};
	auto mat = linear::Matrix<T>(n_row, n_col, vec);

	auto expected = std::vector<T>();
	for (linear::natural c = 0; c < n_col; ++c){
		for (linear::natural r = 0; r < n_row; ++r){
			expected.push_back(vec[c + r*n_col]); 
		}	
	}
		
	size_t i = 0;
	auto assert_increment {true};
	for (auto it = mat.col_begin(); it != mat.col_end(); it++, ++i){
		assert_increment &= *it == expected[i];
	}

	auto assert_decrement {true};
	for (auto it = (--mat.col_end()); it != mat.col_begin(); it--, --i){
		assert_decrement &= *(it) == expected[i-1];
	}

	BOOST_TEST(assert_increment);
}

BOOST_AUTO_TEST_SUITE_END()

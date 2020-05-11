#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <Linear.hpp>
#include <random>
#include <iostream>

using test_types = boost::mpl::list<int, double, float, char>;
using numeric_test_types = boost::mpl::list<double>;

/// Basic tests

BOOST_AUTO_TEST_SUITE(matrix_basic)

BOOST_AUTO_TEST_CASE(matrix_ostream,
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

BOOST_AUTO_TEST_CASE(diag_ctor,
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

BOOST_AUTO_TEST_CASE(vector_ctor,
                      * boost::unit_test::description("Tests matrix constructor from vector of elements") )
{
	
	std::string expected = "1 2 3 4\n5 6 7 8";
	std::vector<int> v {1, 2, 3, 4, 5, 6, 7, 8};

	auto m {linear::Matrix<int>(2, 4, v)};

	std::ostringstream stream;
	stream << m;
	std::string str = stream.str();

	auto check_1 = (str.length() == expected.length());
  	auto check_2 = (std::equal(str.begin(), str.end(), expected.begin()));

	BOOST_TEST(check_1 & check_2);
}

BOOST_AUTO_TEST_CASE(matrix_transpose,
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

BOOST_AUTO_TEST_CASE(inplace_matrix_multiplication,
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

BOOST_AUTO_TEST_CASE(mat_mult_2_2_2_4,
					* boost::unit_test::description("Tests matrix matrix multiplication."))
{
	std::vector<int> m_rhs {1, 2, 3, 4};
	std::vector<int> m_lhs {1, 2, 3, 4, 5, 6, 7, 8}; 
	auto rhs {linear::Matrix<int>(2, 2, m_rhs)};
	auto lhs {linear::Matrix<int>(2, 4, m_lhs)};

	//	1 2	 *	1  2  3  4		=	11	14	 17	 20
	//	4 5		5  6  7  8			23 	30	 37	 44
	//	
	
	auto res = rhs * lhs;
	std::vector<int> expected {11, 14, 17, 20, 23, 30, 37, 44}; 

	int i {0};
	
	auto is_correct {true};
	for (auto it = res.begin(); it != res.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	
	
	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(mat_mult_2_3_3_4,
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

BOOST_AUTO_TEST_CASE(mat_mult_3_3_3_4,
					* boost::unit_test::description("Tests matrix matrix multiplication."))
{
	std::vector<int> m_rhs {1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<int> m_lhs {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}; 
	auto rhs {linear::Matrix<int>(3, 3, m_rhs)};
	auto lhs {linear::Matrix<int>(3, 4, m_lhs)};

	//	1 2 3	 *	1  2  3  4		=	38	44	 50	 56
	//	4 5	6		5  6  7  8			83 	98	113	128
	//	7 8	9		9 10 11 12			128	152	176	200
	
	auto res = rhs * lhs;
	std::vector<int> expected {38, 44, 50, 56, 83, 98, 113, 128, 128, 152, 176, 200}; 

	int i {0};
	
	auto is_correct {true};
	for (auto it = res.begin(); it != res.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	
	
	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(mat_mult_3_3_3_3,
					* boost::unit_test::description("Tests matrix matrix multiplication."))
{
	std::vector<int> m_rhs {1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<int> m_lhs {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
	auto rhs {linear::Matrix<int>(3, 3, m_rhs)};
	auto lhs {linear::Matrix<int>(3, 3, m_lhs)};

	//	1 2 3	 *	1  2  3		=	30	36 	42
	//	4 5	6		4  5  6			66	81	96
	//	7 8	9		7  8  9			102 126 150
	
	auto res = rhs * lhs;
	std::vector<int> expected {30, 36, 42, 66, 81, 96, 102, 126, 150}; 

	int i {0};
	
	auto is_correct {true};
	for (auto it = res.begin(); it != res.end(); ++it, ++i){
		is_correct &= *it == expected[i];
	}	
	
	BOOST_TEST(is_correct);
}

BOOST_AUTO_TEST_CASE(large_mat_mult_1)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(1); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(1.0, 10.0);

	int l = 21;

	using Matrix = linear::Matrix<double>;
	std::vector<double> rhs_m {};
	std::vector<double> lhs_m {};

	for(int i=0; i<l*l; ++i){
		rhs_m.push_back(dis(gen));
		lhs_m.push_back(dis(gen));
	}

	auto rhs {Matrix(l, l, rhs_m)};
	auto lhs {Matrix(l, l, lhs_m)};

	auto res = lhs*rhs;

	//std::cout<< rhs << std::endl;
	//std::cout<< "*" << std::endl;
	//std::cout<< lhs << std::endl;
	//std::cout<< "=" << std::endl;
	//std::cout<< res << std::endl;

	// calculate expected result assuming matrix
	// is row major 
	double max_err = 0.0;
	std::vector<double> expected;
	for(int r=0; r<l; r++){
		for(int c=0; c<l; c++){
			int v = 0;
			for(int k=0; k<l; k++){
				v += lhs_m[r*l + k] * rhs_m[k*l + c];
			}
			auto err = abs(v - res(r, c));
			max_err =  err > max_err? err : max_err;
		}
	}

	BOOST_TEST(max_err < 1.0e-15);
}

BOOST_AUTO_TEST_CASE(large_mat_mult_2)
{
	int l = 9;

	using Matrix = linear::Matrix<double>;
	std::vector<double> rhs_m {};
	std::vector<double> lhs_m {};

	for(int i=0; i<l*l; ++i){
		rhs_m.push_back(1.0);
		lhs_m.push_back(1.0);
	}

	auto rhs {Matrix(l, l, rhs_m)};
	auto lhs {Matrix(l, l, lhs_m)};

	auto res = lhs*rhs;

	//std::cout<< rhs << std::endl;
	//std::cout<< "*" << std::endl;
	//std::cout<< lhs << std::endl;
	//std::cout<< "=" << std::endl;
	//std::cout<< res << std::endl;

	// calculate expected result assuming matrix
	// is row major 
	double max_err = 0.0;
	std::vector<double> expected;
	for(int r=0; r<l; r++){
		for(int c=0; c<l; c++){
			int v = 0;
			for(int k=0; k<l; k++){
				v += lhs_m[r*l + k] * rhs_m[k*l + c];
			}
			auto err = abs(v - res(r, c));
			max_err =  err > max_err? err : max_err;
		}
	}

	BOOST_TEST(max_err < 1.0e-15);
}

BOOST_AUTO_TEST_CASE(mat_inplace_sum,
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

BOOST_AUTO_TEST_CASE(mat_inplace_subtr,
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

BOOST_AUTO_TEST_CASE(mat_instream_exception,
					* boost::unit_test::description("Tests matrix instream operator exception."))
{
	std::vector<int> vec {1, 2, 3, 4, 5, 6, 7};
	linear::Matrix<int> m(3, 2);

	BOOST_CHECK_THROW(m >> vec, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(invalid_vector_ctor,
					 * boost::unit_test::description("Tests matrix constructor from invalid vector of elements")){

	std::vector<int> v {1, 2, 3, 4, 5, 6, 7, 8};
	BOOST_CHECK_THROW(linear::Matrix<int>(2, 5, v), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

/// Templates tests

BOOST_AUTO_TEST_SUITE(matrix_templ)

BOOST_AUTO_TEST_CASE_TEMPLATE(base_ctor, T, test_types){
	
	auto mat = linear::Matrix<T>();
	auto value = linear::natural();

	auto check_1 = mat.nrow() == value;
	auto check_2 = mat.ncol() == value;

	BOOST_TEST(check_1 & check_2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(copy_ctor, T, test_types){
	
	linear::natural n {5};
	auto mat {linear::Matrix<T>(n, n)};
	auto mat_copy {mat};

	BOOST_TEST(not (&mat(0,0) == &mat_copy(0,0)));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(move_op, T, test_types){
	
	linear::natural n {5};
	auto mat {linear::Matrix<T>(n, n)};
	auto ptr = &(mat(0,0));
	
	auto mat_copy = std::move(mat);

	BOOST_TEST(&(mat_copy(0,0)) == ptr);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(row_it, T, test_types){
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

BOOST_AUTO_TEST_CASE_TEMPLATE(row_it_post_incr, T, test_types){
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

BOOST_AUTO_TEST_CASE_TEMPLATE(col_it, T, test_types){

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

BOOST_AUTO_TEST_CASE_TEMPLATE(col_iterator_post_incr, T, test_types){

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

// Numeric templates  

BOOST_AUTO_TEST_SUITE(matrix_num_templ)

BOOST_AUTO_TEST_CASE_TEMPLATE(t_mat_mult, T, numeric_test_types)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    auto seed{rd()};
	std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 10.0);
	std::uniform_int_distribution<> int_dis(1, 50);

	auto lhs_n_row = int_dis(gen);
	auto lhs_n_col = int_dis(gen);
	auto rhs_n_col = int_dis(gen);

	using Matrix = linear::Matrix<T>;
	std::vector<T> rhs_m {};
	std::vector<T> lhs_m {};

	for(int i=0; i<lhs_n_row*lhs_n_col; ++i){
		lhs_m.push_back(static_cast<T>(dis(gen)));
	}

	for(int i=0; i<lhs_n_col*rhs_n_col; ++i){
		rhs_m.push_back(static_cast<T>(dis(gen)));
	}

	auto lhs {Matrix(lhs_n_row, lhs_n_col, lhs_m)};
	auto rhs {Matrix(lhs_n_col, rhs_n_col, rhs_m)};

	auto res = lhs*rhs;
	
	//std::cout<< lhs_n_row << " " << lhs_n_col << std::endl;
	//std::cout<< "*" << std::endl;
	//std::cout<< lhs_n_col << " " << rhs_n_col << std::endl;	
	//std::cout<< "=" << std::endl;
	//std::cout<< res << std::endl;

	// calculate expected result assuming matrix
	// is row major 
	std::vector<T> expected;

	T max_err = 0;
	for(int r=0; r<lhs_n_row; r++){
		for(int c=0; c<rhs_n_col; c++){
			T v = 0;
			for(int k=0; k<lhs_n_col; k++){
				v += lhs_m[r*lhs_n_col + k] * rhs_m[k*rhs_n_col + c];
			}
			auto err = abs(v - res(r, c));
			max_err =  err > max_err? err : max_err;
		}
	}
	
	auto pass{max_err < 1.0e-10};
	if(!pass){
		std::cout << seed << std::endl;
	}

	BOOST_TEST(pass);
}

BOOST_AUTO_TEST_SUITE_END()
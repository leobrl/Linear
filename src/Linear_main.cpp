#include "Linear.hpp"
#include <chrono>
#include <random>

int main(){

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

	using Matrix = linear::Matrix<double>;
	std::vector<double> m {}; 
	int p = 20;
	for(int i=0; i<std::pow(2.0, p);++i){
		m.push_back(dis(gen));
	}

	auto sz = static_cast<linear::natural>(std::pow(2.0, p/2));
	auto rhs {Matrix(sz, sz, m)};
	auto lhs {Matrix(sz, sz, m)};

	auto start = std::chrono::steady_clock::now();
	const int NUMBER_OF_LOOPS{10};
	for(size_t i = 0; i<NUMBER_OF_LOOPS; i++){
		auto temp = rhs * lhs;
	}
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> delta_t = (end-start);

	double avg_time = (delta_t.count())/NUMBER_OF_LOOPS;
	double gflop = 2*(sz*sz*sz)*1.0E-9;

	std::cout << sz << "x"<< sz << ":\n";
	std::cout << "\tAvg time " <<  avg_time << "\n";
	std::cout << "\tGFLOP " << gflop << "\n";
	std::cout << "\tGFLOP/sec " << gflop/avg_time << "\n";
	std::cout << std::endl;

	return 0;
}

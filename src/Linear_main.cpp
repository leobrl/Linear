#include "Linear.hpp"
#include <chrono>
#include <random>


int main(){

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(1, 10);

	std::vector<int> m {}; 
	int p = 20;
	for(int i=0; i<std::pow(2.0, p);++i){
		m.push_back(dis(gen));
	}

	auto sz = static_cast<linear::natural>(std::pow(2.0, p/2));
	auto rhs {linear::Matrix<int>(sz, sz, m)};
	auto lhs {linear::Matrix<int>(sz, sz, m)};

	auto start = std::chrono::steady_clock::now();
	auto res = rhs * lhs;
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;

	std::cout << diff.count() << std::endl;

	return 0;
}

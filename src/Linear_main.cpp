#include "Linear.hpp"
#include <chrono>
#include <random>

int main(){

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(1, 10);

	std::vector<int> m {}; 
	for(int i=0; i<10000;++i){
		m.push_back(dis(gen));
	}

	auto rhs {linear::Matrix<int>(100, 100, m)};
	auto lhs {linear::Matrix<int>(100, 100, m)};


	auto start = std::chrono::steady_clock::now();
	rhs *= lhs;
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;

	std::cout << diff.count() << std::endl;

	return 0;
}

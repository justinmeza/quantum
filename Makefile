all: *.cpp
	clang++ -std=c++20 *.cpp && ./a.out

all: Bra.cpp Bra.hpp Complex.cpp Complex.hpp Const.hpp Const.cpp Ket.cpp Ket.hpp Operator.hpp Operator.cpp quantum.cpp
	clang++ -std=c++20 $^ && ./a.out

test: Bra.cpp Bra.hpp Complex.cpp Complex.hpp Const.hpp Const.cpp Ket.cpp Ket.hpp Operator.hpp Operator.cpp test.cpp
	clang++ -std=c++20 $^ && ./a.out

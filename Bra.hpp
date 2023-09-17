#ifndef __BRA_HPP
#define __BRA_HPP

#include <vector>
#include <iostream>

#include "Complex.hpp"
/* #include "Ket.hpp" */

using namespace std;

class Ket;

class Bra {
	public:
		vector<Complex> elements;

		Bra();
		explicit Bra(vector<Complex> elements);
		/* Bra(const Bra& b); */
		bool operator==(Bra rhs);

		friend Bra operator*(Bra lhs, const double rhs);
		friend Bra operator*(double lhs, Bra rhs);
		friend Bra operator*(Bra lhs, const Complex rhs);
		friend Bra operator*(Complex lhs, Bra rhs);
		friend Complex operator*(Bra lhs, const Ket rhs);

		friend Bra operator+(Bra lhs, const Bra rhs);

		Ket toKet();

	private:
};

ostream& operator<<(ostream& os, const Bra& obj);

#endif // __BRA_HPP

#ifndef __BRA_HPP
#define __BRA_HPP

#include <vector>
#include <iostream>
#include "Complex.hpp"

using namespace std;

class Ket;

class Bra {
	public:
		vector<Complex> elements;

		Bra() {}

		Bra(auto elements) {
			this->elements = elements;
		}

		friend Bra operator*(Bra lhs, const double rhs) {
			Complex scalar(rhs, 0);
			for (auto e = lhs.elements.begin(); e != lhs.elements.end(); e++) {
				*e = *e * scalar;
			}
			return lhs;
		}

		friend Bra operator*(double lhs, Bra rhs) {
			Complex scalar(lhs, 0);
			for (auto e = rhs.elements.begin(); e != rhs.elements.end(); e++) {
				*e = scalar * *e;
			}
			return rhs;
		}

		friend Bra operator*(Bra lhs, const Complex rhs) {
			for (auto e = lhs.elements.begin(); e != lhs.elements.end(); e++) {
				*e = *e * rhs;
			}
			return lhs;
		}

		friend Bra operator*(Complex lhs, Bra rhs) {
			for (auto e = rhs.elements.begin(); e != rhs.elements.end(); e++) {
				*e = lhs * *e;
			}
			return rhs;
		}

		friend Bra operator+(Bra lhs, const Bra rhs) {
			auto el = lhs.elements;
			auto er = rhs.elements;
			if (el.size() != er.size()) {
				cerr << "Attempting to add Bras of different sizes (" << el.size() << " and " << er.size() << ")" << endl;
				exit(1);
			}
			for (auto e = 0; e < el.size(); e++) {
				lhs.elements[e] = el[e] + er[e];
			}
			return lhs;
		}

		friend Complex operator*(Bra lhs, const Ket rhs);

		Ket toKet();

	private:
};

ostream& operator<<(ostream& os, const Bra& obj);

#endif // __BRA_HPP

#ifndef __KET_HPP
#define __KET_HPP

#include <vector>
#include <iostream>

#include "Complex.hpp"

using namespace std;

class Bra;

class Ket {
	public:
		vector<Complex> elements;

		Ket();
		explicit Ket(vector<Complex> elements);
		Ket(const Ket& k);

		bool operator==(Ket rhs) {
			auto eps = 1e-10;
			auto n = 0;
			for (auto e : elements) {
				/* cout << e << " vs. " << rhs.elements[n] << endl; */
				if (e != rhs.elements[n])
					return false;
				n++;
			}
			return true;
		}

		friend Ket operator*(double lhs, Ket rhs) {
			Complex scalar(lhs, 0);
			for (auto e = rhs.elements.begin(); e != rhs.elements.end(); e++) {
				*e = scalar * *e;
			}
			return rhs;
		}

		friend Ket operator*(Ket lhs, const Complex rhs) {
			for (auto e = lhs.elements.begin(); e != lhs.elements.end(); e++) {
				*e = *e * rhs;
			}
			return lhs;
		}

		friend Ket operator*(Complex lhs, Ket rhs) {
			for (auto e = rhs.elements.begin(); e != rhs.elements.end(); e++) {
				*e = lhs * *e;
			}
			return rhs;
		}

		friend Ket operator+(Ket lhs, const Ket rhs) {
			auto el = lhs.elements;
			auto er = rhs.elements;
			if (el.size() != er.size()) {
				cerr << "Attempting to add Kets of different sizes (" << el.size() << " and " << er.size() << ")" << endl;
				exit(1);
			}
			for (auto e = 0; e < el.size(); e++) {
				lhs.elements[e] = el[e] + er[e];
			}
			return lhs;
		}

		Bra toBra();

	private:
};

ostream& operator<<(ostream& os, const Ket& obj);

#endif // __KET_HPP

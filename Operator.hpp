#ifndef __OPERATOR_HPP
#define __OPERATOR_HPP

#include "Complex.hpp"
#include "Bra.hpp"
#include "Ket.hpp"

using namespace std;

class Operator {
	public:
		vector<vector<Complex> > elements;

		Operator() {}

		Operator(auto elements) {
			this->elements = elements;
		}

		friend Ket operator*(Operator lhs, const Ket rhs) {
			Ket ret;
			for (auto r = lhs.elements.begin(); r != lhs.elements.end(); r++) {
				Bra tmp(*r);
				ret.elements.push_back(tmp * rhs);
			}
			return ret;
		}

		friend Bra operator*(Bra lhs, const Operator rhs) {
			Bra ret;
			for (auto c = 0; c < rhs.elements.size(); c++) {
				Ket tmp;
				for (auto r = rhs.elements.begin(); r != rhs.elements.end(); r++) {
					tmp.elements.push_back((*r)[c]);
				}
				ret.elements.push_back(lhs * tmp);
			}
			return ret;
		}

		Operator transpose() {
			Operator ret;
			ret.elements = elements;
			for (auto r = 0; r < elements.size(); r++) {
				for (auto c = 0; c < elements.size(); c++) {
					ret.elements[c][r] = elements[r][c];
				}
			}
			return ret;
		}

		Operator conjugate() {
			Operator ret;
			ret.elements = elements;
			for (auto r = 0; r < elements.size(); r++) {
				for (auto c = 0; c < elements.size(); c++) {
					ret.elements[r][c] = elements[r][c].conjugate();
				}
			}
			return ret;
		}

		Operator hermetian() {
			Operator ret;
			ret.elements = elements;
			return ret.transpose().conjugate();
		}

	private:
};

ostream& operator<<(ostream& os, const Operator& obj)
{
	os << "[";
	for (auto r = 0; r < obj.elements.size(); r++) {
		Bra tmp(obj.elements[r]);
		os << tmp;
		if (r < obj.elements.size() - 1) {
			os << ", ";
		}
	}
	os << "]";
	return os;
}

#endif // __OPERATOR_HPP

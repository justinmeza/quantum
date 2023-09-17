#include "Operator.hpp"

Operator operator*(Ket lhs, Bra rhs) {
	vector<vector<Complex>> ret;
	ret.resize(lhs.elements.size(), vector<Complex>(rhs.elements.size(), 0.0));
	for (auto l = 0; l < lhs.elements.size(); l++) {
		for (auto r = 0; r < rhs.elements.size(); r++) {
			ret[l][r] = lhs.elements[l] * rhs.elements[r];
		}
	}
	return Operator(ret);
}

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

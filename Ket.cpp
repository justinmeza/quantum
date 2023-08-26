#include "Bra.hpp"
#include "Ket.hpp"

Bra Ket::toBra() {
	vector<Complex> ret;
	for (auto e = elements.begin(); e != elements.end(); e++) {
		ret.push_back((*e).conjugate());
	}
	return Bra(ret);
}

ostream& operator<<(ostream& os, const Ket& obj)
{
	os << "[";
	for (auto e = 0; e < obj.elements.size(); e++) {
		os << obj.elements[e];
		if (e != obj.elements.size() - 1) {
			os << ", ";
		}
	}
	os << "]^T";
	return os;
}

#include "Bra.hpp"
#include "Ket.hpp"

Ket Bra::toKet() {
	vector<Complex> ret;
	for (auto e = elements.begin(); e != elements.end(); e++) {
		ret.push_back((*e).conjugate());
	}
	return Ket(ret);
}

Complex operator*(Bra lhs, const Ket rhs) {
	auto el = lhs.elements;
	auto er = rhs.elements;
	Complex acc(0.0, 0.0);
	if (el.size() != er.size()) {
		cerr << "Attempting to multiply Bra and Ket of different sizes (" << el.size() << " and " << er.size() << ")" << endl;
		exit(1);
	}
	for (auto e = 0; e < el.size(); e++) {
		acc = acc + (el[e] * er[e]);
	}
	return acc;
}

ostream& operator<<(ostream& os, const Bra& obj)
{
	os << "[";
	for (auto e = 0; e < obj.elements.size(); e++) {
		os << obj.elements[e];
		if (e != obj.elements.size() - 1) {
			os << ", ";
		}
	}
	os << "]";
	return os;
}

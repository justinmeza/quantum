#include "Bra.hpp"
#include "Complex.hpp"
#include "Ket.hpp"

Bra::Bra() {}

Bra::Bra(vector<Complex> elements) {
	this->elements = elements;
}

/* Bra::Bra(const Bra& b) { */
/* 	this->elements = b.elements; */
/* } */

bool Bra::operator==(Bra rhs) {
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

Bra operator*(Bra lhs, const double rhs) {
	Complex scalar(rhs, 0);
	for (auto e = lhs.elements.begin(); e != lhs.elements.end(); e++) {
		*e = *e * scalar;
	}
	return lhs;
}

Bra operator*(Bra lhs, const Complex rhs) {
	for (auto e = lhs.elements.begin(); e != lhs.elements.end(); e++) {
		*e = *e * rhs;
	}
	return lhs;
}

Bra operator*(Complex lhs, Bra rhs) {
	for (auto e = rhs.elements.begin(); e != rhs.elements.end(); e++) {
		*e = lhs * *e;
	}
	return rhs;
}

Bra operator*(double lhs, Bra rhs) {
	Complex scalar(lhs, 0);
	for (auto e = rhs.elements.begin(); e != rhs.elements.end(); e++) {
		*e = scalar * *e;
	}
	return rhs;
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

Bra operator+(Bra lhs, const Bra rhs) {
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

Ket Bra::toKet() {
	vector<Complex> ret;
	for (auto e = elements.begin(); e != elements.end(); e++) {
		ret.push_back((*e).conjugate());
	}
	return Ket(ret);
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

#include "Complex.hpp"

ostream& operator<<(ostream& os, const Complex& obj)
{
	string sign = obj.getImaginary() > 0.0 ? "+" : "-";
	os << obj.getReal() << " " << sign << " " << abs(obj.getImaginary()) << "i";
	return os;
}

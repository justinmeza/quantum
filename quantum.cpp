#include <vector>
#include <iostream>

#include "Complex.hpp"
#include "Bra.hpp"
#include "Ket.hpp"
#include "Operator.hpp"

using namespace std;

int main(int argc, char** argv) {
	Complex cc(1, 2);
	Complex zero(0, 0);
	Complex one(1, 0);
	Complex two(2, 0);

	cout << "(" << cc << ") + (" << cc << ") = " << cc + cc << endl;
	cout << "(" << cc << ") * (" << cc << ") = " << cc * cc << endl;
	cout << "-3 * (" << cc << ") = " << -3 * cc << endl;
	cout << "(" << cc << ") * 2 = " << cc * 2 << endl;
	cout << "(" << cc << ") - (" << cc << ") = " << cc - cc << endl;
	cout << "(" << cc << ") / (" << cc << ") = " << cc / cc << endl;
	cout << "(" << cc << ") / (" << cc << ")^2 = " << cc / (cc * cc) << endl;
	cout << "(" << cc << ")^* = " << cc.conjugate() << endl;
	cout << "(" << cc << ")^* * (" << cc << ") = " << cc.conjugate() * cc << endl;
	cout << "r^2 = " << cc.getRadius() * cc.getRadius() << endl;

	Ket k(vector<Complex>{cc, cc});
	cout << k << endl;
	cout << k.toBra() << endl;
	cout << k * 2 << endl;
	cout << 3 * k << endl;
	cout << k + k << endl;
	cout << cc * k << endl;
	cout << k * cc << endl;

	Bra b(vector<Complex>{cc, cc});
	cout << b << endl;
	cout << b.toKet() << endl;
	cout << b * 2 << endl;
	cout << 3 * b << endl;
	cout << b + b << endl;
	cout << cc * b << endl;
	cout << b * cc << endl;
	cout << b * cc.conjugate() << endl;
	cout << cc * b.toKet() << endl;

	cout << b * k << endl;
	cout << (k.toBra() * b.toKet()).conjugate() << endl;

	cout << b * (k + k) << endl;
	cout << b * k + b * k << endl;

	Operator M(vector<vector<Complex> >{
		{one, two},
		{two, one}
	});
	cout << M << endl;
	cout << M * k << endl;
	cout << M * (k + k) << endl;
	cout << M * k + M * k << endl;
	cout << M * (2 * k) << endl;
	cout << b * M << endl;
	cout << M.transpose() << endl;
	cout << M.conjugate() << endl;
	cout << "Hermetian = " << M.hermetian() << endl;
	cout << M * k << endl;
	cout << k.toBra() * M.hermetian() << endl;
	cout << (M * k).toBra() << endl;

	return 0;
}

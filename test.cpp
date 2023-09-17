#include <iostream>

#include "Complex.hpp"
#include "Const.hpp"
#include "Bra.hpp"
#include "Ket.hpp"
#include "Operator.hpp"

static int n = 0;
static int f = 0;

#define TEST(s, x) { \
	bool p = true; \
	IS_TRUE(x); \
	cout << "[" << (p ? "PASSED" : "FAILED") << "] " << n << ": " << s << endl; \
	n++; \
}

#define IS_TRUE(x) { \
	if (!(x)) { \
		cout << __FUNCTION__ << " failed on line " << __LINE__ << endl; \
		p = false; \
		f++; \
	} \
}

using namespace std;

void test() {
	Complex cc(1, 2);
	/* Complex zero(0, 0); */
	/* Complex one(1, 0); */
	/* Complex two(2, 0); */
	/* Complex minus_one(-1, 0); */
	/* Complex i(0, 1); */
	/* Complex minus_i(0, -1); */
	/* Complex oneoverroottwo(1 / sqrt(2), 0); */

	TEST("complex/onetwo/real",
			cc.getReal() == 1);
	TEST("complex/onetwo/imaginary",
			cc.getImaginary() == 2);
	TEST("complex/zero/real",
			Const::zero.getReal() == 0);
	TEST("complex/zero/imaginary",
			Const::zero.getImaginary() == 0);
	TEST("complex/one/real",
			Const::one.getReal() == 1);
	TEST("complex/one/imaginary",
			Const::one.getImaginary() == 0);
	TEST("complex/two/real",
			(Const::one * 2).getReal() == 2);
	TEST("complex/two/imaginary",
			(Const::one * 2).getImaginary() == 0);
	TEST("complex/minus_one/real",
			Const::minus_one.getReal() == -1);
	TEST("complex/minus_one/imaginary",
			Const::minus_one.getImaginary() == 0);
	TEST("complex/i/real",
			Const::i.getReal() == 0);
	TEST("complex/i/imaginary",
			Const::i.getImaginary() == 1);
	TEST("complex/minus_i/real",
			Const::minus_i.getReal() == 0);
	TEST("complex/minus_i/imaginary",
			Const::minus_i.getImaginary() == -1);
	TEST("complex/root_two_inv/real",
			Const::root_two_inv.getReal() == 1 / sqrt(2));
	TEST("complex/root_two_inv/imaginary",
			Const::root_two_inv.getImaginary() == 0);

	/* cout << "(" << cc << ") + (" << cc << ") = " << cc + cc << endl; */
	TEST("complex/add",
			(cc + cc) == Complex(2, 4));
	/* cout << "(" << cc << ") * (" << cc << ") = " << cc * cc << endl; */
	TEST("complex/multiply",
			(cc * cc) == Complex(-3, 4));
	/* cout << "-3 * (" << cc << ") = " << -3 * cc << endl; */
	TEST("complex/scalar_prefix",
			(-3 * cc) == Complex(-3, -6));
	/* cout << "(" << cc << ") * 2 = " << cc * 2 << endl; */
	TEST("complex/scalar_postfix",
			(cc * 2) == Complex(2, 4));
	/* cout << "(" << cc << ") - (" << cc << ") = " << cc - cc << endl; */
	TEST("complex/subtract",
			(cc - cc) == Complex(0, 0));
	/* cout << "(" << cc << ") / (" << cc << ") = " << cc / cc << endl; */
	TEST("complex/divide",
			(cc / cc) == Complex(1, 0));
	/* cout << "(" << cc << ") / (" << cc << ")^2 = " << cc / (cc * cc) << endl; */
	TEST("complex/divide_fraction",
			(cc / (cc * cc)) == Complex(0.2, -0.4));
	/* cout << "(" << cc << ")^* = " << cc.conjugate() << endl; */
	TEST("complex/conjugate",
			cc.conjugate() == Complex(1, -2));
	/* cout << "(" << cc << ")^* * (" << cc << ") = " << cc.conjugate() * cc << endl; */
	TEST("complex/conjugate_times",
			(cc.conjugate() * cc) == Complex(5, 0));

	Ket k(vector<Complex>{cc, cc});
	Bra b(vector<Complex>{cc, cc});

	/* cout << k << endl; */
	TEST("ket/elements",
			(k.elements[0] == cc)
			&& (k.elements[1] == cc));
	/* cout << k.toBra() << endl; */
	TEST("ket/to_bra",
			(k.toBra().elements[0] == cc.conjugate())
			&& (k.toBra().elements[1] == cc.conjugate()));
	/* cout << k * 2 << endl; */
	TEST("ket/scalar_postfix",
			(k * 2).elements[0] == (cc * 2)
			&& (k * 2).elements[1] == (cc * 2));
	/* cout << 3 * k << endl; */
	TEST("ket/scalar_prefix",
			(3 * k).elements[0] == (3 * cc)
			&& (3 * k).elements[1] == (3 * cc));
	/* cout << k + k << endl; */
	TEST("ket/addition",
			((k + k).elements[0] == (cc + cc))
			&& ((k + k).elements[1] == (cc + cc)));
	/* cout << cc * k << endl; */
	TEST("ket/multiply_complex_prefix",
			((cc * k).elements[0] == (cc * cc))
			&& ((cc * k).elements[1] == (cc * cc)));
	/* cout << k * cc << endl; */
	TEST("ket/multiply_complex_postfix",
			((k * cc).elements[0] == (cc * cc))
			&& ((k * cc).elements[1] == (cc * cc)));

	/* cout << b << endl; */
	TEST("bra/elements",
			(b.elements[0] == cc)
			&& (b.elements[1] == cc));
	/* cout << b.toKet() << endl; */
	TEST("bra/to_ket",
			(b.toKet().elements[0] == cc.conjugate())
			&& (b.toKet().elements[1] == cc.conjugate()));
	/* cout << b * 2 << endl; */
	TEST("bra/scalar_postfix",
			(b * 2).elements[0] == (cc * 2)
			&& (b * 2).elements[1] == (cc * 2));
	/* cout << 3 * b << endl; */
	TEST("bra/scalar_prefix",
			(3 * b).elements[0] == (3 * cc)
			&& (3 * b).elements[1] == (3 * cc));
	/* cout << b + b << endl; */
	TEST("bra/addition",
			((b + b).elements[0] == (cc + cc))
			&& ((b + b).elements[1] == (cc + cc)));
	/* cout << cc * b << endl; */
	TEST("bra/multiply_complex_prefix",
			((cc * b).elements[0] == (cc * cc))
			&& ((cc * b).elements[1] == (cc * cc)));
	/* cout << b * cc << endl; */
	TEST("bra/multiply_complex_postfix",
			((b * cc).elements[0] == (cc * cc))
			&& ((b * cc).elements[1] == (cc * cc)));
	/* cout << b * cc.conjugate() << endl; */
	/* cout << cc * b.toKet() << endl; */

	/* cout << b * k << endl; */
	TEST("braket/multiply",
			((b * k) == (cc * cc + cc * cc)));
	/* cout << (k.toBra() * b.toKet()).conjugate() << endl; */
	TEST("braket/from_conversion",
			((k.toBra() * b.toKet()).conjugate() == (cc * cc + cc * cc)));

	/* cout << b * (k + k) << endl; */
	TEST("braket/distribution",
			((b * (k + k)) == (b * k + b * k)));
	/* cout << b * k + b * k << endl; */

	Operator M(vector<vector<Complex> >{
		{Const::one,     Const::i},
		{Const::one * 4, Const::one * 3}
	});
	/* cout << M << endl; */
	TEST("operator/definition",
			((M.elements[0][0] == Const::one)
			 && (M.elements[1][1] == Const::one * 3)
			 && (M.elements[0][1] == Const::i)
			 && (M.elements[1][0] == Const::one * 4)));
	/* cout << M * k << endl; */
	TEST("operator/ket_multiplication",
			(((M * k).elements[0] == (k.elements[0] * Const::one + k.elements[0] * Const::i))
			 && ((M * k).elements[1] == (k.elements[1] * Const::one * 4 + k.elements[1] * Const::one * 3))));
	/* cout << M * (k + k) << endl; */
	/* cout << M * k + M * k << endl; */
	/* cout << M * (2 * k) << endl; */
	/* cout << b * M << endl; */
	TEST("operator/bra_multiplication",
			(((b * M).elements[0] == (b.elements[0] * Const::one + b.elements[0] * Const::one * 4))
			 && ((b * M).elements[1] == (b.elements[1] * Const::i + b.elements[1] * Const::one * 3))));
	/* cout << M.transpose() << endl; */
	TEST("operator/transpose",
			(M.transpose().elements[0][0] == Const::one
			 && M.transpose().elements[0][1] == Const::one * 4
			 && M.transpose().elements[1][0] == Const::i
			 && M.transpose().elements[1][1] == Const::one * 3));
	/* cout << M.conjugate() << endl; */
	TEST("operator/conjugate",
			(M.conjugate().elements[0][0] == Const::one
			 && M.conjugate().elements[0][1] == -1 * Const::i
			 && M.conjugate().elements[1][0] == Const::one * 4
			 && M.conjugate().elements[1][1] == Const::one * 3));
	cout << "Hermetian = " << M.hermetian() << endl;
	cout << M * k << endl;
	cout << k.toBra() * M.hermetian() << endl;
	cout << (M * k).toBra() << endl;
}

int main(int argc, char** argv) {
	test();
	cout << f << " FAILED / " << (n - f) << " PASSED / " << n << " TOTAL" << endl;
}

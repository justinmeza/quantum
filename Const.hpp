#ifndef __CONST_HPP
#define __CONST_HPP

#include "Complex.hpp"
#include "Bra.hpp"
#include "Ket.hpp"
#include "Operator.hpp"

class Const {
	public:
		static const double h, h_bar;
		static const Complex zero, one, minus_one, i, minus_i, root_two_inv, minus_root_two_inv, i_root_two_inv;
		static const Operator pauli_x, pauli_y , pauli_z;
		static const Ket up, down, left, right, in, out;
};

#endif

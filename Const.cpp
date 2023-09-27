#include "Const.hpp"

const double Const::h = 1.054571726e-34;
const double Const::h_bar = Const::h / (2 * M_PI);

const Complex Const::zero(0, 0);
const Complex Const::one(1, 0);
const Complex Const::minus_one(-1, 0);
const Complex Const::i(0, 1);
const Complex Const::minus_i(0, -1);
const Complex Const::root_two_inv(1 / sqrt(2), 0);
const Complex Const::minus_root_two_inv(-1 / sqrt(2), 0);
const Complex Const::i_root_two_inv(0, 1 / sqrt(2));

const Operator Const::pauli_x(vector<vector<Complex> >{
		{Const::zero, Const::one},
		{Const::one,  Const::zero}});
const Operator Const::pauli_y(vector<vector<Complex> >{
		{Const::zero, Const::minus_i},
		{Const::i,    Const::zero}});
const Operator Const::pauli_z(vector<vector<Complex> >{
		{Const::one,  Const::zero},
		{Const::zero, Const::minus_one}});

const Ket Const::up(vector<Complex>{
		Const::one,
		Const::zero});
const Ket Const::down(vector<Complex>{
		Const::zero,
		Const::one});
const Ket Const::right(vector<Complex>{
		Const::root_two_inv,
		Const::root_two_inv});
const Ket Const::left(vector<Complex>{
		Const::root_two_inv,
		Const::minus_root_two_inv});
const Ket Const::in(vector<Complex>{
		Complex(1.0 / sqrt(2.0), 0),
		Complex(0, 1.0 / sqrt(2.0))});
const Ket Const::out(vector<Complex>{
		Complex(1.0 / sqrt(2.0), 0),
		Complex(0, -1.0 / sqrt(2.0))});

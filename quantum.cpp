#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "Complex.hpp"
#include "Const.hpp"
#include "Bra.hpp"
#include "Ket.hpp"
#include "Operator.hpp"

using namespace std;

// Pattern for building a time-dependent Schrödinger Ket.
void schrodingerKet(Operator hamiltonian, Ket initial_state, double time) {
	double h = 1.054571726e-34;
	double h_bar = h / (2 * M_PI);
	cout << "┅┅┅┅┅ Schrödinger Ket ┅┅┅┅┅" << endl;

	// Define the Hamiltonian operator.
	Operator H = hamiltonian;

	// Define the initial state, |ψ(0)⟩.
	Ket Phi_0 = initial_state;

	// Find the eigenvalues and eigenvectors of H.
	map<Complex, Ket> eigen = H.eigen();

	// Confirm that H|λ_n⟩ = λ_n|λ_n⟩.
	for (auto e : eigen) {
		cout << H * e.second << " ≟ " << e.first * e.second << endl;
	}

	// Calculate the initial coefficients ɑ_n(0).
	vector<Complex> alpha_0;
	for (auto e : eigen) {
		alpha_0.push_back(e.second.toBra() * Phi_0);
	}

	// Derive the time-dependent state |ψ(t)⟩.
	Ket Phi_t;
	int n = 0;
	for (auto e : eigen) {
		// Convert from Euler form to Cartesian form.
		Complex c(
				cos(-1.0 / h_bar * e.first.getReal() * time),
				sin(-1.0 / h_bar * e.first.getReal() * time)
		);
		Ket k = alpha_0[n] * c * e.second;
		if (n == 0)
			Phi_t = k;
		else
			Phi_t = Phi_t + k;
		n++;
	}
	cout << "|ψ(" << time << ")⟩ = " << Phi_t << endl;

	// Compute the probabilities P_λ(t).
	map<double, double> P_lamda;
	for (auto e : eigen) {
		P_lamda.insert(make_pair(
				e.first.getReal(),
				((e.second.toBra() * Phi_t).abs() * (e.second.toBra() * Phi_t).abs()).getReal()
		));
	}
	for (auto p : P_lamda) {
		cout << "P_" << p.first << "(" << time << ") = " << p.second << endl;
	}
}

int main(int argc, char** argv) {

	// Mastering Quantum Mechanics:  Essentialsl, Theory, and Applications
	// Beam splitter 1 (BS1).
	Operator U_1(vector<vector<Complex>>{
		{-1 * Const::root_two_inv, Const::root_two_inv},
		{Const::root_two_inv,      Const::root_two_inv}
	});
	// Beam splitter 2 (BS2).
	Operator U_2(vector<vector<Complex>>{
		{Const::root_two_inv, Const::root_two_inv},
		{Const::root_two_inv, -1 * Const::root_two_inv}
	});
	// Input wave function (incident photon state).
	Ket ab(vector<Complex>{Const::zero, Const::one});
	cout << "---------" << endl;
	cout << U_2 * U_1 << endl;
	// Output wave function (photon detector state).
	Ket D = U_2 * U_1 * ab;
	cout << "D0 = " << D.elements[0] << endl;
	cout << "D1 = " << D.elements[1] << endl;
	// Blocking the lower path.
	Ket BS1_out = U_1 * ab;
	cout << BS1_out << endl;
	Ket ab_blocked = BS1_out;
	ab_blocked.elements[1] = 0;
	cout << ab_blocked << endl;
	Ket BS2_out = U_2 * ab_blocked;
	cout << BS2_out << endl;

	// High reflexivity beam splitter.
	/* double N = 1.0; */
	/* double N = 0.5; */
	double N = 1e12;
	Operator U_hr(vector<vector<Complex>>{
		{Complex(cos(M_PI / (2 * N)), 0.0), Complex(0.0, sin(M_PI / (2 * N)))},
		{Complex(0.0, sin(M_PI / (2 * N))), Complex(cos(M_PI / (2 * N)), 0.0)}
	});
	cout << U_hr * U_hr.hermetian() << endl;
	cout << U_hr * ab << endl;
	// Reflexivity.
	double R = pow(cos(M_PI / (2 * N)), 2.0);
	// Transmissivity.
	double T = pow(sin(M_PI / (2 * N)), 2.0);
	cout << "R = " << R << endl;
	cout << "T = " << T << endl;

	// Plank's constant.
	double h = 1.054571726e-34;
	double h_bar = h / (2 * M_PI);
	double omega = 1.0;  // TODO:  What shoud this be?

	schrodingerKet(
			(omega * h_bar / 2.0) * Const::pauli_z,
			Const::up,
			1.0
	);

	return 0;
}

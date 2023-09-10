#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "Complex.hpp"
#include "Bra.hpp"
#include "Ket.hpp"
#include "Operator.hpp"

using namespace std;

Operator realBlockHermetian(Operator o) {
	int len = o.elements.size();
	Operator ret;
	ret.elements.resize(len * 2, vector<Complex>());
	for (auto r = 0; r < len * 2; r++) {
		ret.elements[r].resize(len * 2, Complex());
		for (auto c = 0; c < len * 2; c++) {
			if (r < len && c < len) {
				// A
				ret.elements[r][c] = o.elements[r][c].getReal();
			} else if (r > len - 1 && c < len) {
				// B
				ret.elements[r][c] = o.elements[r - len][c].getImaginary();
			} else if (r < len && c > len - 1) {
				// -B
				ret.elements[r][c] = -1 * o.elements[r][c - len].getImaginary();
			} else if (r > len - 1 && c > len - 1) {
				// A
				ret.elements[r][c] = o.elements[r - len][c - len].getReal();
			}
		}
	}
	return ret;
}

void jacobiReal(
		const vector<vector<double>>& matrix,
		vector<double>& eigenvalues,
		vector<vector<double>>& eigenvectors,
		int maxIterations = 1000,
		double tolerance = 1e-6
) {
	int n = matrix.size();
	eigenvalues.resize(n);
	eigenvectors.resize(n, vector<double>(n, 0.0));
	vector<vector<double>> A = matrix;

	// Initializing eigenvectors as identity matrix
	for(int i = 0; i < n; ++i) {
		eigenvectors[i][i] = 1.0;
	}

	for(int iter = 0; iter < maxIterations; ++iter) {
		double maxOffDiagonal = 0.0;
		int p, q;
		for(int i = 0; i < n; ++i) {
			for(int j = i+1; j < n; ++j) {
				if(fabs(A[i][j]) > maxOffDiagonal) {
					maxOffDiagonal = fabs(A[i][j]);
					p = i;
					q = j;
				}
			}
		}

		if(maxOffDiagonal < tolerance) {
			for(int i = 0; i < n; ++i) {
				eigenvalues[i] = A[i][i];
			}
			return;
		}

		double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
		double t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
		if(theta < 0.0) {
			t = -t;
		}
		double c = 1.0 / sqrt(1.0 + t*t);
		double s = t * c;
		double tau = s / (1.0 + c);

		double apq = A[p][q];
		A[p][q] = 0.0;
		A[p][p] -= t * apq;
		A[q][q] += t * apq;

		for(int i = 0; i < n; ++i) {
			if(i != p && i != q) {
				double aip = A[i][p];
				double aiq = A[i][q];
				A[i][p] = A[p][i] = aip - s * (aiq + tau * aip);
				A[i][q] = A[q][i] = aiq + s * (aip - tau * aiq);
			}
		}

		for(int i = 0; i < n; ++i) {
			double uip = eigenvectors[i][p];
			double uiq = eigenvectors[i][q];
			eigenvectors[i][p] = uip - s * (uiq + tau * uip);
			eigenvectors[i][q] = uiq + s * (uip - tau * uiq);
		}
	}

	// If reached here, maximum iterations exceeded
	cerr << "Jacobi method did not converge within the maximum number of iterations" << endl;
}

void jacobi(
		int n,
		vector<vector<Complex> > &a,
		vector<vector<Complex> > &v,
		vector<Complex> &d,
		int nrot,
		double EPS
) {
	int i, j, ip, iq;
	Complex tresh, theta, tau, t, sm, s, h, g, c;
	vector<Complex> b, z; b.resize(n); z.resize(n);
	for (ip = 0; ip < n; ip++) {
		for (iq = 0; iq < n; iq++) {
			v[ip][iq] = 0.0;
			v[ip][ip] = 1.0;
		}
	}
	for (ip = 0; ip < n; ip++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	for (i = 1; i <= 50; i++) {
		sm = 0.0;
		for (ip = 0; ip < n - 1; ip++) {
			for (iq = ip + 1; iq < n; iq++) {
				sm = sm + a[ip][iq].abs();
			}
		}
		// NOTE:  Correct assumption?  .abs()?
		if (sm.getReal() == 0.0)
			return;
		if (i < 4)
			tresh = 0.2 * sm / (n * n);
		else
			tresh = 0.0;
		for (ip = 0; ip < n - 1; ip++) {
			for (iq = ip + 1; iq < n; iq++) {
				g = 100.0 * a[ip][iq].abs();
				// NOTE:  Correct assumption?  .abs()?
				if (i > 4 && g.getReal() <= (EPS * d[ip].abs()).getReal() && g.getReal() <= (EPS * d[iq].abs()).getReal())
					a[ip][iq] = 0.0;
				// NOTE:  Correct assumption?  .abs()?
				else if (a[ip][iq].abs().getReal() > tresh.getReal()) {
					h = d[iq] - d[ip];
					if (g.getReal() <= (EPS * h.abs()).getReal())
						t = a[ip][iq] / h;
					else {
						theta = 0.5 * h / a[ip][iq];
						t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
						// NOTE:  Correct assumption?  .abs()?
						if (theta.getReal() < 0.0) t = -1 * t;
					}
					c = 1.0 / (1 + t * t).sqrt();
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];
					z[ip] = z[ip] - h;
					z[iq] = z[iq] + h;
					d[ip] = d[ip] - h;
					d[iq] = d[iq] + h;
					a[ip][iq] = 0.0;
					for (j = 0; j < ip; j++) {
						Complex x = a[j][ip];
						Complex y = a[j][iq];
						a[j][ip] = x - s * (y + x * tau);
						a[j][iq] = y + s * (x - y * tau);
					}
					for (j = ip + 1; j < iq; j++) {
						Complex x = a[ip][j];
						Complex y = a[j][iq];
						a[ip][j] = x - s * (y + x * tau);
						a[j][iq] = y + s * (x - y * tau);
					}
					for (j = iq + 1; j < n; j++) {
						Complex x = a[ip][j];
						Complex y = a[iq][j];
						a[ip][j] = x - s * (y + x * tau);
						a[iq][j] = y + s * (x - y * tau);
					}
					for (j = 0; j < n; j++) {
						Complex x = v[j][ip];
						Complex y = v[j][iq];
						v[j][ip] = x - s * (y + x * tau);
						v[j][iq] = y + s * (x - y * tau);
					}
					++nrot;
				}
			}
		}
		for (ip = 0; ip < n; ip++) {
			b[ip] = b[ip] + z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	cerr << "Jacobi method did not converge within the maximum number of iterations." << endl;
}

int main(int argc, char** argv) {
	Complex cc(1, 2);
	Complex zero(0, 0);
	Complex one(1, 0);
	Complex two(2, 0);
	Complex minus_one(-1, 0);
	Complex i(0, 1);
	Complex minus_i(0, -1);
	Complex oneoverroottwo(1/sqrt(2), 0);

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

	// Pauli matrices
	Operator sigma_z(vector<vector<Complex> >{
		{one, zero},
		{zero, minus_one}
	});
	Operator sigma_x(vector<vector<Complex> >{
		{zero, one},
		{one, zero}
	});
	Operator sigma_y(vector<vector<Complex> >{
		{zero, minus_i},       // Row-major.
		{i, zero}
	});
	cout << "σ_z = " << sigma_z << endl;
	cout << "σ_x = " << sigma_x << endl;
	cout << "σ_y = " << sigma_y << endl;

	Ket up(vector<Complex>{one, zero});
	Ket down(vector<Complex>{zero, one});

	Ket right(vector<Complex>{Complex(1.0 / sqrt(2.0), 0), Complex(1.0 / sqrt(2.0), 0)});
	Ket left(vector<Complex>{Complex(1.0 / sqrt(2.0), 0), Complex(-1.0 / sqrt(2.0), 0)});

	Ket in(vector<Complex>{Complex(1.0 / sqrt(2.0), 0), Complex(0, 1.0 / sqrt(2.0))});
	Ket out(vector<Complex>{Complex(1.0 / sqrt(2.0), 0), Complex(0, -1.0 / sqrt(2.0))});

	cout << sigma_z * up << endl;
	cout << up << endl;
	cout << sigma_z * down << endl;
	cout << -1.0 * down << endl;
	cout << sigma_z * right << endl;
	cout << right.toBra() * sigma_z * right << endl;

	vector<vector<double>> matrix = {
		{2.0, -1.0, 0.0},
		{-1.0, 2.0, -1.0},
		{0.0, -1.0, 2.0}
	};

	Operator sigma_y_real(vector<vector<Complex> >{
		{zero, zero, zero, one},       // Row-major.
		{zero, zero, minus_one, zero},
		{zero, minus_one, zero, zero},
		{one, zero, zero, zero}
	});

	Operator sigma_y_block = realBlockHermetian(sigma_y);
	cout << "Real block Hermetian: " << sigma_y_block << endl;

	vector<Complex> eigenvalues;
	vector<vector<Complex>> eigenvectors;

	map<Complex, Ket> eigen_y = sigma_y.eigen();
	for (auto e : eigen_y) {
		cout << sigma_y * e.second * -1 << " =?= " << e.first * e.second << endl;
	}

	/* sigma_y_block.jacobi(sigma_y_block.elements, eigenvalues, eigenvectors); */
	/* cout << "Eigenvalues: "; */
	/* for (const auto& val : eigenvalues) { */
	/* 	cout << val << " "; */
	/* } */
	/* cout << endl; */
	/* cout << "Eigenvectors: " << endl; */
	/* for (const auto& row : eigenvectors) { */
	/* 	for (const auto& val : row) { */
	/* 		cout << val << " "; */
	/* 	} */
	/* 	cout << endl; */
	/* } */

	/* Ket jin(eigenvectors[0]); */
	/* Ket jout(eigenvectors[1]); */

	/* cout << jin * jin.toBra() << endl; */
	/* cout << jout * jout.toBra() << endl; */

	// Mastering Quantum Mechanics:  Essentialsl, Theory, and Applications
	// Beam splitter 1 (BS1).
	Operator U_1(vector<vector<Complex>>{
		{-1 * oneoverroottwo, oneoverroottwo},
		{oneoverroottwo,      oneoverroottwo}
	});
	// Beam splitter 2 (BS2).
	Operator U_2(vector<vector<Complex>>{
		{oneoverroottwo,      oneoverroottwo},
		{oneoverroottwo, -1 * oneoverroottwo}
	});
	// Input wave function (incident photon state).
	Ket ab(vector<Complex>{zero, one});
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

	/* jacobi(2, sigma_y.elements, eigenvectors, eigenvalues, 0, 1e-16); */
	/* jacobiMethodComplex(sigma_y_real.elements, eigenvalues, eigenvectors); */
	/* cout << "Eigenvalues: "; */
	/* for (const auto& val : eigenvalues) { */
	/* 	cout << val << " "; */
	/* } */
	/* cout << endl; */

	/* cout << "Eigenvectors: " << endl; */
	/* for (const auto& row : eigenvectors) { */
	/* 	for (const auto& val : row) { */
	/* 		cout << val << " "; */
	/* 	} */
	/* 	cout << endl; */
	/* } */

	return 0;
}

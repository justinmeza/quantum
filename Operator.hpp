#ifndef __OPERATOR_HPP
#define __OPERATOR_HPP

#include <map>

#include "Complex.hpp"
#include "Bra.hpp"
#include "Ket.hpp"

using namespace std;

class Operator {
	public:
		vector<vector<Complex>> elements;  // Row-major.

		Operator() {}

		explicit Operator(auto elements) {
			for (auto e : elements) {
				this->elements.push_back(e);
			}
		}

		friend Operator operator*(Complex lhs, const Operator rhs) {
			Operator ret;
			for (auto r : rhs.elements) {
				Bra b;
				for (auto c : r) {
					b.elements.push_back(c);
				}
				ret.elements.push_back(b.elements);
			}
			return rhs;
		}

		friend Ket operator*(Operator lhs, const Ket rhs) {
			Ket ret;
			for (auto r = lhs.elements.begin(); r != lhs.elements.end(); r++) {
				Bra tmp(*r);
				ret.elements.push_back(tmp * rhs);
			}
			return ret;
		}

		friend Bra operator*(Bra lhs, const Operator rhs) {
			Bra ret;
			for (auto c = 0; c < rhs.elements.size(); c++) {
				Ket tmp;
				for (auto r = rhs.elements.begin(); r != rhs.elements.end(); r++) {
					tmp.elements.push_back((*r)[c]);
				}
				ret.elements.push_back(lhs * tmp);
			}
			return ret;
		}

		friend Operator operator*(Operator lhs, const Operator rhs) {
			int rows = lhs.elements.size();
			int cols = rhs.elements[0].size();
			Operator ret;
			ret.elements.resize(rows, vector<Complex>(cols, Complex(0.0, 0.0)));
			for (auto r = 0; r < rows; r++) {
				for (auto c = 0; c < cols; c++) {
					Ket col;
					col.elements.resize(cols, Complex(0.0, 0.0));
					for (auto x = 0; x < cols; x++) {
						col.elements[x] = rhs.elements[x][c];
					}
					Bra row;
					row.elements = lhs.elements[r];
					ret.elements[r][c] = row * col;
				}
			}
			return ret;
		}

		Operator transpose() {
			Operator ret;
			ret.elements = elements;
			for (auto r = 0; r < elements.size(); r++) {
				for (auto c = 0; c < elements.size(); c++) {
					ret.elements[c][r] = elements[r][c];
				}
			}
			return ret;
		}

		Operator conjugate() {
			Operator ret;
			ret.elements = elements;
			for (auto r = 0; r < elements.size(); r++) {
				for (auto c = 0; c < elements.size(); c++) {
					ret.elements[r][c] = elements[r][c].conjugate();
				}
			}
			return ret;
		}

		Operator hermetian() {
			Operator ret;
			ret.elements = elements;
			return ret.transpose().conjugate();
		}

		Operator realBlock() const {
			int len = elements.size();
			Operator ret;
			ret.elements.resize(len * 2, vector<Complex>());
			for (auto r = 0; r < len * 2; r++) {
				ret.elements[r].resize(len * 2, Complex());
				for (auto c = 0; c < len * 2; c++) {
					if (r < len && c < len) {
						// A
						ret.elements[r][c] = elements[r][c].getReal();
					} else if (r > len - 1 && c < len) {
						// B
						ret.elements[r][c] = elements[r - len][c].getImaginary();
					} else if (r < len && c > len - 1) {
						// -B
						ret.elements[r][c] = -1 * elements[r][c - len].getImaginary();
					} else if (r > len - 1 && c > len - 1) {
						// A
						ret.elements[r][c] = elements[r - len][c - len].getReal();
					}
				}
			}
			return ret;
		}

		map<Complex, Ket> eigen() const {
			Operator o = this->realBlock();
			vector<Complex> eigenvalues;
			vector<vector<Complex>> eigenvectors;
			this->jacobi(o.elements, eigenvalues, eigenvectors);
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
			// Match eigenvalues to correct eigenvectors.
			map<Complex, Ket> ret;
			for (auto n = 0; n < eigenvalues.size(); n++) {
				Complex eigval = eigenvalues[n];
				for (auto m : eigenvectors) {
					Ket u, v;
					for (auto l = 0; l < m.size() / 2; l++) {
						u.elements.push_back(m[l]);
						v.elements.push_back(m[l + (m.size() / 2)]);
					}
					Ket eigvec(u + Complex(0.0, 1.0) * v);
					if ((((*this) * eigvec) == (eigval * eigvec)) && (ret.find(eigval) == ret.end())) {
						/* cout << eigval << ": " << eigvec << endl; */
						ret.insert(pair<Complex, Ket>(eigval, eigvec));
						break;
					}
				}
			}
			return ret;
		}

		void jacobi(
				const vector<vector<Complex>>& A,
				vector<Complex>& eigenvalues,
				vector<vector<Complex>>& eigenvectors,
				int maxIterations = 1000,
				double tolerance = 1e-6,
				bool debug = false
		) const {
			int n = A.size();
			eigenvalues.resize(n);
			eigenvectors.resize(n, vector<Complex>(n, Complex(0.0, 0.0)));

			// Initialize eigenvectors as identity matrix.
			for (int i = 0; i < n; ++i) {
				eigenvectors[i][i] = Complex(1.0, 0.0);
			}

			vector<vector<Complex>> B = A;

			for (int iter = 0; iter < maxIterations; ++iter) {
				// Find the largest off-diagonal element and its indices.
				double maxOffDiagonal = 0.0;
				int p, q;
				for (int i = 0; i < n; ++i) {
					for (int j = i + 1; j < n; ++j) {
						if (B[i][j].abs().getReal() > maxOffDiagonal) {
							maxOffDiagonal = B[i][j].abs().getReal();
							p = i;
							q = j;
						}
					}
				}
				if (debug) {
					cout << "Max off-diagonal: " << maxOffDiagonal << endl;
					cout << "p = " << p << ", " << "q = " << q << endl;
				}

				/* cout << "Iteration " << iter << ": " << Operator(B) << endl; */
				// Check for convergence.
				if (maxOffDiagonal < tolerance) {
					for (int i = 0; i < n; ++i) {
						eigenvalues[i] = B[i][i];
					}
					if (debug) {
						cout << "Jacobi method converged in " << iter << " iterations." << endl;
					}
					return;
				}

				// Compute the Jacobi rotation parameters.
				Complex theta = (B[q][q] - B[p][p]) / (2.0 * B[p][q]);
				Complex t = 1.0 / (theta.abs() + (1.0 + theta.abs() * theta.abs()).sqrt());
				if (theta.getReal() < 0.0) {
					t = -1 * t;
				}
				Complex c = 1.0 / (1.0 + t.abs() * t.abs()).sqrt();
				Complex s = t * c;

				// Apply the Jacobi rotation to B and the eigenvectors matrix.
				Complex temp;
				Complex conjugate_s = s.conjugate();
				for (int j = 0; j < n; ++j) {
					if (j != p && j != q) {
						temp = c * B[p][j] - s * B[q][j];
						B[q][j] = s * B[p][j] + c * B[q][j];
						B[p][j] = temp;

						// Also update the symmetric elements.
						B[j][p] = temp.conjugate();
						B[j][q] = B[q][j].conjugate();
					}

					// Update the eigenvector matrix.
					temp = c * eigenvectors[j][p] - conjugate_s * eigenvectors[j][q];
					eigenvectors[j][q] = s * eigenvectors[j][p] + c * eigenvectors[j][q];
					eigenvectors[j][p] = temp;
				}

				// Update the diagonal and off-diagonal elements.
				Complex App = c * c * B[p][p] - 2.0 * c * s * B[p][q] + s * s * B[q][q];
				Complex Aqq = s * s * B[p][p] + 2.0 * c * conjugate_s * B[p][q] + c * c * B[q][q];
				B[p][q] = B[q][p] = (c * c - s * s) * B[p][q] + c * s * (B[p][p] - B[q][q]);
				B[p][p] = App;
				B[q][q] = Aqq;
			}

			// If reached here, maximum iterations exceeded
			cerr << "Jacobi method did not converge within the maximum number of iterations" << endl;
		}

	private:
};

Operator operator*(Ket lhs, Bra rhs);
ostream& operator<<(ostream& os, const Operator& obj);

#endif // __OPERATOR_HPP

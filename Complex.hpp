#ifndef __COMPLEX_HPP
#define __COMPLEX_HPP

#include <vector>
#include <iostream>

using namespace std;

class Complex {
	public:
		Complex(double x, double y) {
			this->x = x;
			this->y = y;
		}

		double getReal() const {
			return this->x;
		}

		void setReal(double x) {
			this->x = x;
		}

		double getImaginary() const {
			return this->y;
		}

		void setImaginary(double y) {
			this->y = y;
		}

		friend Complex operator+(Complex lhs, const Complex& rhs) {
			lhs.setReal(lhs.getReal() + rhs.getReal());
			lhs.setImaginary(lhs.getImaginary() + rhs.getImaginary());
			return lhs;
		}

		friend Complex operator+(Complex lhs, const double rhs) {
			Complex scalar(rhs, 0);
			return lhs + scalar;
		}

		friend Complex operator+(double lhs, const Complex rhs) {
			Complex scalar(lhs, 0);
			return scalar + rhs;
		}

		friend Complex operator*(Complex lhs, const double rhs) {
			Complex scalar(rhs, 0);
			return lhs * scalar;
		}

		friend Complex operator*(double lhs, const Complex rhs) {
			Complex scalar(lhs, 0);
			return scalar * rhs;
		}

		friend Complex operator*(Complex lhs, const Complex& rhs) {
			double x = lhs.getReal() * rhs.getReal() - lhs.getImaginary() * rhs.getImaginary();
			double y = lhs.getReal() * rhs.getImaginary() + lhs.getImaginary() * rhs.getReal();
			lhs.setReal(x);
			lhs.setImaginary(y);
			return lhs;
		}

		friend Complex operator-(Complex lhs, const Complex& rhs) {
			return lhs + (-1 * rhs);
		}

		friend Complex operator-(double lhs, const Complex rhs) {
			Complex scalar(lhs, 0);
			return scalar - rhs;
		}

		friend Complex operator-(Complex lhs, const double rhs) {
			Complex scalar(rhs, 0);
			return lhs - scalar;
		}

		friend Complex operator/(Complex lhs, const Complex& rhs) {
			double x = (lhs.getReal() * rhs.getReal() + lhs.getImaginary() * rhs.getImaginary())
					/ (rhs.getReal() * rhs.getReal() + rhs.getImaginary() * rhs.getImaginary());
			double y = (lhs.getImaginary() * rhs.getReal() - lhs.getReal() * rhs.getImaginary())
					/ (rhs.getReal() * rhs.getReal() + rhs.getImaginary() * rhs.getImaginary());
			lhs.setReal(x);
			lhs.setImaginary(y);
			return lhs;
		}

		friend Complex operator/(Complex lhs, const double rhs) {
			Complex scalar(rhs, 0);
			return lhs / scalar;
		}

		friend Complex operator/(double lhs, const Complex rhs) {
			Complex scalar(lhs, 0);
			return scalar / rhs;
		}

		Complex conjugate() {
			return Complex(this->x, -1 * this->y);
		}

		double getRadius() {
			return sqrt(x * x + y * y);
		}

		double getTheta() {
			return atan(y / x);
		}

	private:
		double x;
		double y;
};

ostream& operator<<(ostream& os, const Complex& obj);

#endif // __COMPLEX_HPP

#ifndef __COMPLEX_HPP
#define __COMPLEX_HPP

#include <vector>
#include <iostream>

using namespace std;

class Complex {
	public:
		Complex(double x = 0.0, double y = 0.0) {
			this->x = x;
			this->y = y;
		}

		bool operator <(const Complex& rhs) const
    {
        return (this->getReal() < rhs.getReal() && this->getImaginary() < rhs.getImaginary())
					|| (this->getReal() == rhs.getReal() && this->getImaginary() < rhs.getImaginary())
					|| (this->getReal() < rhs.getReal() && this->getImaginary() == rhs.getImaginary());
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
			double x = lhs.getReal() + rhs.getReal();
			double y = lhs.getImaginary() + rhs.getImaginary();
			return Complex(x, y);
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
			return Complex(x, y);
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
			return Complex(x, y);
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
			return ::sqrt(x * x + y * y);
		}

		double getTheta() {
			return atan(y / x);
		}

		Complex abs() const {
			return Complex(::sqrt(x * x + y * y), 0.0);
		}

		Complex sqrt() {
			/* int sign = y < 0 ? -1 : 1; */
			/* return Complex(::sqrt(0.5 * (::sqrt(x * x + y * y) + x)), sign * ::sqrt(0.5 * (::sqrt(x * x + y * y) - x))); */
			double r = this->abs().getReal();
			return ::sqrt(r) * (*this + r) / (*this + r).abs();

			/* double r = ::sqrt(::hypot(this->getReal(), this->getImaginary())); */
			/* double theta = ::atan2(this->getImaginary(), this->getReal()); */
			/* double sqrt_r = ::sqrt(r); */
			/* double angle = theta / 2; */
			/* double real = sqrt_r * ::cos(angle); */
			/* double imag = sqrt_r * ::sin(angle); */
			/* return Complex(real, imag); */
		}

	private:
		double x;
		double y;
};

ostream& operator<<(ostream& os, const Complex& obj);

#endif // __COMPLEX_HPP

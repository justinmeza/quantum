# Quantum

A simple, quantum mechanics sandbox.  I created this as a learning aid while
reading Susskind and Friedman's _Quantum Mechanics: The Theoretical Minimum_
and Zwiebath's _Mastering Quantum Mechanics: Essentials, Theory, and
Applications_.

The main goal was to have a standalone implementation of the basic mathematical
concepts involved in solving quantum mechanical systems.  While I rely on
`<cmath>` for some basic trigonometric functionality, I wrote my own simple
complex number implementation so that the reader can understand how things like
the norm or root of a complex number is calculated.

Another design decision was to avoid off-the-shelf Eigensystem solvers in favor
of using the basic approach described in Sections 11.5 in _Numerical Recipes:
The Art of Scientific Computing_, with a bit of help from GPT-4 to streamline
some of the code.  I wanted to avoid copyright stuff and use precisely what I
needed to solve the problem at hand (an eigensystem for small Hermetian
matrices).

## Example: Time-Dependent Schrödinger Solution

Here is an example of using this library to solve the time-dependent
Schrödinger equation for the most basic quantum system consisting of a single
spin starting in a particular state from Section 4.13 in _Theoretical Minimum_:

```
void schrodingerKet(Operator hamiltonian, Ket initial_state, double time) {
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
                cos(-1.0 / Const::h_bar * e.first.getReal() * time),
                sin(-1.0 / Const::h_bar * e.first.getReal() * time)
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
```

Then, we can use the function like so:

```
double omega = 1.0;
double time = 1.0;
Ket initial = Const::up;

schrodingerKet(
    (omega * Const::h_bar / 2.0) * Const::pauli_z,
    initial,
    time
);
```

Which outputs what we expect:  If the initial state of the quantum system is with an `up` spin and we measure it along the z-axis, the probability of us measuring `up` (corresponding to the Eigenvalue `1`) is `1.0`.

```
P_-1(1) = 0
P_1(1) = 1
```

A more interesting example is if we change the initial state to be `right` and then perform the same measurement along the z-axis:

```
P_-1(1) = 0.5
P_1(1) = 0.5
```

Since our direction of measurement is orthogonal to the initial state, we are equally likely to measure either `up` or `down`.

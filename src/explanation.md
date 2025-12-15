# Explaining the Lane-Emden Equation Solver

This document provides an explanation of the Python script that solves the Lane-Emden equation using the Adomian Decomposition Method (ADM). This explanation is intended for an audience with a mathematical background.

## The Lane-Emden Equation

The Lane-Emden equation is a second-order ordinary differential equation that describes the structure of a self-gravitating, spherically symmetric, polytropic fluid. It is a fundamental equation in astrophysics, used to model stars and other celestial objects. The equation is given by:

$$
\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2\frac{d\theta}{d\xi}\right) + \theta^n = 0
$$

with initial conditions:

$$ 
\theta(0) = 1, \quad \theta'(0) = 0
$$ 

where $\xi$ is a dimensionless radius and $\theta$ is related to the density (and thus pressure) by $\rho = \rho_c \theta^n$, where $\rho_c$ is the central density. The parameter $n$ is the polytropic index.

## The Adomian Decomposition Method (ADM)

The Adomian Decomposition Method (ADM) is a semi-analytical method for solving ordinary and partial nonlinear differential equations. The method was developed by George Adomian. The key idea of ADM is to decompose the solution into a series of functions and to decompose the nonlinear term into a series of polynomials called Adomian polynomials.

Consider a general nonlinear ordinary differential equation in the form:

$$Lu + Ru + Nu = g$$

where $L$ is the highest order derivative which is assumed to be invertible, $R$ is the remainder of the linear operator, $Nu$ represents the nonlinear term, and $g$ is the source term.

Applying the inverse operator $L^{-1}$ to the equation, we get:

$$u = \phi + L^{-1}(g) - L^{-1}(Ru) - L^{-1}(Nu)$$

where $\phi$ is the solution to $Lu=0$ that satisfies the given initial conditions.

The ADM assumes the solution $u$ can be decomposed into an infinite series:

$$u = \sum_{n=0}^{\infty} u_n$$

and the nonlinear term $Nu$ can be decomposed into an infinite series of Adomian polynomials:

$$Nu = \sum_{n=0}^{\infty} A_n$$

where $A_n$ are the Adomian polynomials of $u_0, u_1, \dots, u_n$. The Adomian polynomials are calculated by the formula:

$$A_n = \frac{1}{n!}\frac{d^n}{d\lambda^n}\left[N\left(\sum_{i=0}^{\infty} \lambda^i u_i\right)\right]_{\lambda=0}$$

Substituting the series decompositions into the equation for $u$, we get:

$$\sum_{n=0}^{\infty} u_n = \phi + L^{-1}(g) - L^{-1}\left(R\left(\sum_{n=0}^{\infty} u_n\right)\right) - L^{-1}\left(\sum_{n=0}^{\infty} A_n\right)$$

From this, we can identify the components of the solution:

$$u_0 = \phi + L^{-1}(g)$$

$$u_{n+1} = -L^{-1}(Ru_n) - L^{-1}(A_n)$$

This provides a recursive scheme to calculate the terms of the solution series.

## Implementation in Python

The provided Python script implements the ADM to solve a specific form of the Lane-Emden equation. The script uses the `sympy` library for symbolic mathematics, which allows for exact and analytical computations.

### The `LaneEmdenSolver` class

The `LaneEmdenSolver` class encapsulates the logic for solving the Lane-Emden equation.

#### Initialization

The constructor `__init__` takes the number of terms `num_terms` to calculate for the series solution.

#### The Inverse Operator

The `_inverse_operator` method implements the inverse operator $L^{-1}$. For the Lane-Emden equation, the operator $L$ is $L = \frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2\frac{d}{d\xi}\right)$. The inverse operator $L^{-1}$ is a double integral operator:

$$L^{-1}(f) = \int_0^\xi s^{-2} \int_0^s t^2 f(t) dt ds$$

#### Adomian Polynomials

The `_calculate_adomian_polynomial` method calculates the Adomian polynomials for the nonlinear term. The script considers a specific form of the Lane-Emden equation where the nonlinear term is $N(\theta) = (\theta^2 - C)^{1.5}$.

#### The `solve` method

The `solve` method orchestrates the solution process. It starts with the initial solution component $\theta_0 = 1$ (from the initial conditions). Then, it iteratively calculates the next terms of the series using the recursive formula:

$$\theta_{k+1} = -L^{-1}(A_k)$$

where $A_k$ is the $k$-th Adomian polynomial. This provides a recursive scheme to calculate the terms of the solution series.

### The `main` function

The `main` function creates an instance of the `LaneEmdenSolver` and calls the `solve` method to compute the solution. The script is set to calculate 7 terms of the series solution.

## Results

The script prints the first few terms of the series solution for $\theta(\xi)$. The output of the script is the symbolic expression for the approximate solution. For example, the first few terms are:

$$\theta_0=1$$
$$\theta_1=-\frac{\xi^2}{6}(1-C)^{3/2}$$
$$\dots$$

The final approximate solution is the sum of these terms:

$$\theta(\xi) = \sum_{i=0}^{N-1} \theta_i$$

where $N$ is the number of terms.

This provides a powerful way to obtain an analytical approximation of the solution to the Lane-Emden equation, which is particularly useful for cases where an exact analytical solution is not known.

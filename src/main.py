#!/usr/bin/env python

"""
# SPDX-FileType: SOURCE
# SPDX-FileType: APPLICATION
# SPDX-FileCopyrightText: © 2025 Matihus Alberth Molina Larios <matihus.molina@unmsm.edu.pe> Carlos Alonso Aznarán Laos <caznaranl@uni.pe>
# SPDX-License-Identifier: GPL-2.0-or-later
This script calculates the approximate solution of the Lane-Emden equation
using the Adomian Decomposition Method (ADM).
Based on the work "Adomian Decomposition Method" by Juan Pablo Ballén López and Pablo Alberto León Velasco.
Retrieved from https://repository.eafit.edu.co/server/api/core/bitstreams/7c769b30-4d00-4f73-985d-df1a256e5fa8/content
"""

from sympy import simplify
from sympy.abc import c, lamda, s, t, xi
from sympy.core import Integer
from sympy.core import diff as D
from sympy.functions.combinatorial.factorials import factorial
from sympy.integrals import integrate
from sympy.printing import pprint


def solve_lane_emden_adomian(num_terms):
    """
    Calculates the approximate solution of the Lane-Emden equation using ADM.

    Parameters:
    - num_terms (int): The number of series terms to calculate.

    Returns:
    - approx_solution: The sum of the calculated terms.
    - components: A list containing each individual theta_k.
    """

    print(f"--- Solving Lane-Emden for (y^2-C)^1.5 with {num_terms} terms ---\n")

    # List to store the components theta_0, theta_1, ...
    thetas = []

    # Initial step: theta_0 = 1 (based on the problem's initial conditions)
    theta_0 = Integer(1)
    thetas.append(theta_0)
    print(f"theta_0 = {theta_0}")

    # Helper function for the inverse operator L^(-1)
    # L^(-1)(f) = Integral_0^x ( s^(-2) * Integral_0^s ( t^2 * f(t) ) dt ) ds
    def inverse_operator(func, variable):
        # Inner integral: int(t^2 * f(t), t, 0, s)
        inner_integral = integrate(t**2 * func.subs(variable, t), (t, 0, s))
        # Outer integral: int(s^(-2) * inner_integral, s, 0, x)
        return integrate(s ** (-2) * inner_integral, (s, 0, variable))

    # Iteration to find theta_1, theta_2...
    for k in range(num_terms - 1):
        # A. Construct the Adomian Polynomial A_k for f(u) = u^m
        # Definition: A_k = (1/k!) * d^k/d(lamda)^k [ (sum(u_i*lamda^i))^m ] | lamda=0

        # Build the parametrized partial sum u(lamda)
        parametrized_sum = sum(thetas[i] * lamda**i for i in range(len(thetas)))

        # Apply the non-linearity (u^2 - c)^1.5
        non_linear_function = (parametrized_sum**2 - c) ** 1.5

        # Differentiate k times with respect to lamda
        derivative = D(non_linear_function, lamda, k)

        # Evaluate at lamda = 0 and divide by k!
        A_k = derivative.subs(lamda, 0) / factorial(k)

        # B. Calculate the next term: theta_{k+1} = - L^(-1)( A_k )
        next_theta = -inverse_operator(A_k, xi)
        next_theta = simplify(next_theta)

        thetas.append(next_theta)
        print(f"theta_{k + 1} = {next_theta}")

    # Construct the final solution
    approx_solution = sum(thetas)
    print("\nApproximate solution (first terms):")
    pprint(approx_solution)

    return approx_solution, thetas


sol_gen, _ = solve_lane_emden_adomian(7)

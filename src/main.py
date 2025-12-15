#!/usr/bin/env python
"""
# SPDX-FileType: SOURCE
# SPDX-FileType: APPLICATION
# SPDX-FileCopyrightText: © 2025 Matihus Alberth Molina Larios <matihus.molina@unmsm.edu.pe>, Carlos Alonso Aznarán Laos <caznaranl@uni.pe>
# SPDX-License-Identifier: GPL-2.0-or-later
This script calculates the approximate solution of the Lane-Emden equation
using the Adomian Decomposition Method (ADM).
Based on the work "Adomian Decomposition Method" by Juan Pablo Ballén López and Pablo Alberto León Velasco.
Retrieved from https://repository.eafit.edu.co/server/api/core/bitstreams/7c769b30-4d00-4f73-985d-df1a256e5fa8/content
"""

from typing import List, Tuple

from sympy.abc import lamda, m, q, s, t, u, xi
from sympy.core import Expr, Integer, Rational, Symbol
from sympy.core import diff as D
from sympy.core.numbers import pi
from sympy.functions.combinatorial.factorials import factorial
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.hyperbolic import cosh, sinh
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.integrals import integrate
from sympy.printing import pprint, pretty
from sympy.simplify import simplify


class LaneEmdenSolver:
    """Solves the Lane-Emden equation using the Adomian Decomposition Method (ADM).

    This class encapsulates the state and steps of the ADM to find an
    approximate series solution for a specific form of the Lane-Emden equation.
    """

    def __init__(
        self, num_terms: int, non_linear_func: Expr, initial_condition: int = 1
    ):
        """
        Initializes the solver.

        Args:
            num_terms: The number of terms to calculate for the series solution.
            non_linear_func: The non-linear function f(u) in the equation, where u is a Symbol('u').
            initial_condition: The initial condition u(0) for the equation.
        """
        if not isinstance(num_terms, int) or num_terms <= 0:
            raise ValueError("Number of terms must be a positive integer.")
        self.num_terms = num_terms
        self.non_linear_func = non_linear_func
        self.initial_condition = Integer(initial_condition)
        self.u = u
        self.us: List[Expr] = []
        self.u_symbols: List[Symbol] = []  # Added to store symbols for printing

    def _inverse_operator(self, func: Expr, variable: Expr) -> Expr:
        """
        Calculates the inverse operator L^(-1).

        L^(-1)(f) = Integral_0^x ( s^(-2) * Integral_0^s ( t^2 * f(t) ) dt ) ds

        Args:
            func: The SymPy expression to apply the operator to.
            variable: The independent variable of the function (e.g., xi).

        Returns:
            The result of the inverse operation as a SymPy expression.
        """
        inner_integral = integrate(t**2 * func.subs(variable, t), (t, 0, s))
        return integrate(s ** (-2) * inner_integral, (s, 0, variable))

    def _calculate_adomian_polynomial(self, k: int) -> Expr:
        """
        Constructs the Adomian Polynomial A_k for the non-linearity.

        The non-linearity is provided by the user.
        The polynomial is defined as:
        A_k = (1/k!) * d^k/d(lamda)^k [ f(sum(s_i*lamda^i)) ] | lamda=0

        Args:
            k: The index of the Adomian polynomial to calculate.

        Returns:
            The Adomian polynomial A_k as a SymPy expression.
        """
        parametrized_sum = sum(self.us[i] * lamda**i for i in range(len(self.us)))
        non_linear_function = self.non_linear_func.subs(self.u, parametrized_sum)
        derivative = D(non_linear_function, lamda, k)
        return derivative.subs(lamda, 0) / factorial(k)

    def solve(self) -> Tuple[Expr, List[Expr]]:
        """
        Calculates the approximate solution of the Lane-Emden equation.

        Returns:
            A tuple containing:
            - The symbolic expression for the approximate solution.
            - A list of the individual solution components (us).
        """
        print(
            f"--- Solving Lane-Emden for f(u) = {self.non_linear_func} with {self.num_terms} terms ---\
"
        )
        self.us = []
        self.u_symbols = []

        # Initial step: u_0 = 1 (from initial conditions)
        u_0_sym = Symbol("u_0")
        self.u_symbols.append(u_0_sym)
        u_0_val = self.initial_condition
        self.us.append(u_0_val)
        print(
            f"{pretty(u_0_sym, use_unicode=True)} = {pretty(u_0_val, use_unicode=True)}"
        )

        # Iteratively find the remaining terms
        for k in range(self.num_terms - 1):
            adomian_poly = self._calculate_adomian_polynomial(k)
            next_u_val = -self._inverse_operator(adomian_poly, xi)
            next_u_val = simplify(next_u_val)
            self.us.append(next_u_val)
            next_u_sym = Symbol(f"u_{k + 1}")
            self.u_symbols.append(next_u_sym)
            print(
                f"{pretty(next_u_sym, use_unicode=True)} = {pretty(next_u_val, use_unicode=True)}"
            )

        approx_solution = sum(self.us)
        print("\nApproximate solution (first terms):")
        pprint(approx_solution, use_unicode=True)

        return approx_solution, self.us


def main():
    """
    Main function to run the Lane-Emden solver.
    """
    try:
        LaneEmdenSolver(num_terms=6, non_linear_func=u**m, initial_condition=1).solve()
        LaneEmdenSolver(
            num_terms=6,
            non_linear_func=(u**2 - (q**2 + 1)) ** Rational(3, 2),
            initial_condition=1,
        ).solve()
        LaneEmdenSolver(
            num_terms=6, non_linear_func=exp(u), initial_condition=0
        ).solve()
        LaneEmdenSolver(
            num_terms=6, non_linear_func=exp(-u), initial_condition=0
        ).solve()
        LaneEmdenSolver(
            num_terms=6, non_linear_func=sin(u), initial_condition=1
        ).solve()
        LaneEmdenSolver(
            num_terms=6, non_linear_func=cos(u), initial_condition=1 + pi / 2
        ).solve()
        LaneEmdenSolver(
            num_terms=6, non_linear_func=sinh(u), initial_condition=1
        ).solve()
        LaneEmdenSolver(
            num_terms=6, non_linear_func=cosh(u), initial_condition=1
        ).solve()
    except ValueError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()

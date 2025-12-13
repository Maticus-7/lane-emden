#!/usr/bin/env python
import antigravity
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def resolver_lane_emden_adomian(n_terminos):
    """
    Calcula la solución aproximada de la ecuación de Lane-Emden usando ADM.

    Parámetros:
    - indice_m (int/float/Symbol): El índice politrópico 'm' (puede ser número o símbolo).
    - n_terminos (int): Número de términos de la serie a calcular.

    Retorna:
    - solucion_aprox: La suma de los términos calculados.
    - componentes: Lista con cada theta_k individual.
    """
    # 1. Definir variables
    xi = sp.Symbol('xi')
    lam = sp.Symbol('lambda')
    C = sp.Symbol('C')

    print(f"--- Solucionando Lane-Emden para (y^2-C)^1.5 con {n_terminos} términos ---\n")

    # Lista para guardar las componentes theta_0, theta_1, ...
    thetas = []

    # 2. Paso inicial: theta_0 = 0 (por la condición inicial y(0)=0)
    theta_0 = sp.Integer(1)
    thetas.append(theta_0)
    print(f"theta_0 = {theta_0}")

    # Función auxiliar para el operador inverso L^(-1)
    # L^(-1)(f) = Integral_0^x ( s^(-2) * Integral_0^s ( t^2 * f(t) ) dt ) ds
    def operador_inverso(func, variable):
        t = sp.Symbol('t', dummy=True)
        s = sp.Symbol('s', dummy=True)
        # Integral interna: int(t^2 * f(t), t, 0, s)
        int_interna = sp.integrate(t **2 * func.subs(variable, t), (t, 0, s))
        # Integral externa: int(s^(-2) * int_interna, s, 0, x)
        return sp.integrate(s**(-2) * int_interna, (s, 0, variable))

    # 3. Iteración para encontrar theta_1, theta_2...
    for k in range(n_terminos - 1):
        # A. Construir el Polinomio de Adomian A_k para f(u) = u^m
        # Usamos la definición: A_k = (1/k!) * d^k/d(lam)^k [ (sum(u_i*lam^i))^m ] | lam=0

        # Construir la suma parcial u(lambda)
        suma_parametrizada = sum(thetas[i] * lam**i for i in range(len(thetas)))

        # Aplicar la no-linealidad u^m
        funcion_no_lineal = (suma_parametrizada**2-C)**1.5

        # Derivar k veces respecto a lambda
        derivada = sp.diff(funcion_no_lineal, lam, k)

        # Evaluar en lambda = 0 y dividir por k!
        A_k = derivada.subs(lam, 0) / sp.factorial(k)

        # B. Calcular el siguiente término: theta_{k+1} = - L^(-1)( A_k )
        siguiente_theta = - operador_inverso(A_k, xi)
        siguiente_theta = sp.simplify(siguiente_theta)

        thetas.append(siguiente_theta)
        print(f"theta_{k+1} = {siguiente_theta}")

    # 4. Construir la solución final
    solucion_aprox = sum(thetas)
    print("\nSolución aproximada (primeros términos):")
    sp.pprint(solucion_aprox)

    return solucion_aprox, thetas



# --- EJEMPLOS DE USO ---

xi = sp.Symbol('xi') # Define xi for global use in plot_comparison calls

# Caso 3: e^y simbólico (Solución general)
sol_gen, _ = resolver_lane_emden_adomian(7)


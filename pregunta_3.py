import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def lagrange(x_vals, y_vals):
    x = sp.Symbol('x')
    n = len(x_vals)
    L = 0
    for i in range(n):
        term = y_vals[i]
        for j in range(n):
            if j != i:
                term *= (x - x_vals[j]) / (x_vals[i] - x_vals[j])
        L += term
    return sp.simplify(L)

# Datos de la Pregunta 2
y_vals = [1.0, 0.69696245, -0.00758192, -0.50450963, -0.44956304, 0.0]
x_vals = [1, 2, 3, 4, 5, 6]

# Obtener el polinomio de interpolación
polinomio = lagrange(x_vals, y_vals)
print("Polinomio de Lagrange simplificado:")
sp.pprint(polinomio)

# Convertir a función numérica para graficar
f_lagrange = sp.lambdify(sp.Symbol('x'), polinomio, modules=['numpy'])

# Graficar
x_plot = np.linspace(1, 6, 500)
y_plot = f_lagrange(x_plot)

plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_plot, label="Polinomio de Lagrange", color='blue')
plt.plot(x_vals, y_vals, 'ro', label="Datos originales (Problema 2)")
plt.title("Interpolación de Lagrange vs. Datos Originales")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()
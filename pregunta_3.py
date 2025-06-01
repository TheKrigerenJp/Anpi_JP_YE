import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# ====================================================
# Función para construir el polinomio de Lagrange
# Dado un conjunto de puntos (x_vals, y_vals)
# ====================================================
def lagrange(x_vals, y_vals):
    x = sp.Symbol('x')   # Variable simbólica 
    n = len(x_vals)   # Número de puntos 
    L = 0   # Inicializacion del polinomio 

    # Construcción del polinomio usando la formula de Lagrange
    for i in range(n):
        term = y_vals[i]   # Comenzando con y_i 
        for j in range(n):
            if j != i:
                # Producto del cociente de (x-x_j) / (x_i-x_j)
                term *= (x - x_vals[j]) / (x_vals[i] - x_vals[j])
        L += term   # Suma cada término L_i(x)*y_i

    
    return sp.simplify(L) # Simplificamos el polinomio simbolico final

# ====================================================
# Datos del problema (Pregunta 2)
# ====================================================
y_vals = [1.0, 0.69696245, -0.00758192, -0.50450963, -0.44956304, 0.0]
x_vals = [1, 2, 3, 4, 5, 6]

# ====================================================
# Generación del polinomio de interpolación de Lagrange
# ====================================================
polinomio = lagrange(x_vals, y_vals)

#Imprimimos el polinomio de forma legible
print("Polinomio de Lagrange simplificado:")
sp.pprint(polinomio)

# ====================================================
# Conversión del polinomio simbólico a función numérica
# para evaluación rápida en arreglos de NumPy
# ====================================================
f_lagrange = sp.lambdify(sp.Symbol('x'), polinomio, modules=['numpy'])

# ====================================================
# Graficar la función interpolada
# ====================================================
x_plot = np.linspace(1, 6, 500) # Intervalo de evaluación
y_plot = f_lagrange(x_plot) # Evalución del polinomio 

plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_plot, label="Polinomio de Lagrange", color='blue') # Curva interpolada 
plt.plot(x_vals, y_vals, 'ro', label="Datos originales (Problema 2)") # Puntos originales
plt.title("Interpolación de Lagrange vs. Datos Originales")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()

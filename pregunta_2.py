import numpy as np
import matplotlib.pyplot as plt


# Método de Thomas para resolver sistemas tridiagonales
def metodoDeThomas(a, b, c, d):
    n = len(d)
    # Copias para no modificar los arreglos originales
    a, b, c, d = map(np.array, (a, b, c, d))
    for i in range(1, n):
        w = a[i - 1] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    x = np.zeros(n)
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x


# Función que implementa el método de diferencias finitas para EDOs de segundo orden
def edo2(p, q, r, h, a, b, y0, yn):
    n = int((b - a) / h)
    x = np.linspace(a, b, n + 1)

    A = np.zeros(n - 1)
    B = np.zeros(n - 1)
    C = np.zeros(n - 1)
    D = np.zeros(n - 1)

    for i in range(1, n):
        xi = x[i]
        A[i - 1] = (-h/2)*p(xi)-1
        B[i - 1] = 2+(h**2)*q(xi)
        C[i - 1] = (h/2)*p(xi)-1
        D[i - 1] = (-h**2)*r(xi)

    # Ajustar extremos por condiciones de frontera
    D[0] += ((h/2)*p(x[1])+1) * y0
    D[n-2] += ((-h/2)*p(x[n])+1) * yn

    y_inner = metodoDeThomas(A[1:], B, C[:-1], D)
    y_total = np.concatenate(([y0], y_inner, [yn]))
    print(y_total)
    return x, y_total


# Funciones dadas por el problema
def p(x): return -1 / x


def q(x): return (1 / (4 * x ** 2)) - 1


def r(x): return 0


# Solución exacta del problema
def y_exacta(x):
    return np.sin(6 - x) / (np.sin(5) * np.sqrt(x))


# Valores de h a probar
h_values = [1, 0.5, 0.2, 0.1, 0.01]

# Graficar soluciones
plt.figure(figsize=(10, 6))
for h in h_values:
    x, y_aprox = edo2(p, q, r, h, 1, 6, 1, 0)
    plt.plot(x, y_aprox, label=f'h = {h}')

# Solución exacta
x_exact = np.linspace(1, 6, 500)
plt.plot(x_exact, y_exacta(x_exact), 'k--', label='Solución exacta')

plt.title('Aproximación por Diferencias Finitas vs Solución Exacta')
plt.xlabel('x')
plt.ylabel('y(x)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

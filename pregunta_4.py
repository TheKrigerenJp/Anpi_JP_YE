import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# Función original del problema
# =====================================================
def f(x):
    return np.sin(6 - x) / (np.sin(5) * np.sqrt(x))

# =====================================================
# Datos tabulados del problema
# =====================================================
x_vals = [1, 2, 3, 4, 5, 6]
y_vals = [1.0, 0.69696245, -0.00758192, -0.50450963, -0.44956304, 0.0]

# =====================================================
# Método de Thomas para resolver sistemas tridiagonales
# =====================================================
def metodoDeThomas(a, b, c, d):
    n = len(d)
    # Copias de seguridad para no modificar los originales
    a, b, c, d = map(np.array, (a, b, c, d))

    #Eliminación hacia adelante 
    for i in range(1, n):
        w = a[i - 1] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]

    # Sustitución hacia atrás 
    x = np.zeros(n)
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x

# =====================================================
# Generación de los trazadores cúbicos naturales
# =====================================================
def trazadorCubico(x, y):
    n = len(x)
    y = np.array(y)
    h = np.diff(x) # h_i = x_{i+1} - x_i  

    # Construcción del sistema tridiagonal
    a = h[:-1]   # Subdiagonal
    b = 2 * (h[:-1] + h[1:])   # DIagonal principal 
    c = h[1:]   # Superdiagonal
    d = 6 * ((y[2:] - y[1:-1]) / h[1:] - (y[1:-1] - y[:-2]) / h[:-1])   # Término independiente

    # Resolvemos el sistema usando el metodo de Thomas
    m_inner = metodoDeThomas(a, b, c, d)

    # Incorporar las condiciones naturales (m_0 = m_n = 0)
    m = np.zeros(n)
    m[1:-1] = m_inner

    # Calculamos los coeficientes para cada polinomio cúbico
    splines = []
    for i in range(n - 1):
        hi = h[i]
        A = (m[i + 1] - m[i]) / (6 * hi)
        B = m[i] / 2
        C = (y[i + 1] - y[i]) / hi - hi * (2 * m[i] + m[i + 1]) / 6
        D = y[i]
        splines.append((A, B, C, D, x[i])) # Se almacena los coeficientes junto con x_i
    return splines

# =====================================================
# Evaluación del trazador cúbico en puntos intermedios
# =====================================================
def evaluarTrazador(splines, x_eval):
    y_eval = []
    for x in x_eval:
        for A, B, C, D, xi in splines:
            if xi <= x <= xi + 1:  # If para validar para h = 1
                dx = x - xi
                y_eval.append(A * dx ** 3 + B * dx ** 2 + C * dx + D)
                break
    return np.array(y_eval)

# =====================================================
# Aplicación del trazador cúbico natural y comparación
# =====================================================
splines = trazadorCubico(x_vals, y_vals)
x_fine = np.linspace(x_vals[0], x_vals[-1], 500)
y_spline = evaluarTrazador(splines, x_fine)
y_original = f(x_fine)

# =====================================================
# Gráfica comparativa
# =====================================================
plt.figure(figsize=(10, 6))
plt.plot(x_fine, y_original, label="Función original", color='black')
plt.plot(x_fine, y_spline, label="Trazador cúbico natural", color='green')
plt.plot(x_vals, y_vals, 'ro', label="Puntos de interpolación")
plt.title("Interpolación con Trazadores Cúbicos Naturales")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()

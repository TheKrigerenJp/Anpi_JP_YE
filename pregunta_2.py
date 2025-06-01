import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------
# Método de Thomas:
# Algoritmo para resolver sistemas tridiagonales
# -----------------------------------------------------
def metodoDeThomas(a, b, c, d):
    n = len(d)
    
    # Se convierte las listas a arrays de NumPy para operaciones eficientes
    a, b, c, d = map(np.array, (a, b, c, d))

    # For para hacer la eliminación hacia adelante
    for i in range(1, n):
        # Calculamos el multiplicador que elimina el elemento a[i-1]
        w = a[i - 1] / b[i - 1]

        # Luego modificamos el pivote b[i] y el término independiente d[i]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    
    # Para continuar hacemos ahora la sustitución hacia atrás
    x = np.zeros(n)
    # En la siguiente linea se calcula el valor de la solución
    x[-1] = d[-1] / b[-1]

    # Y el siguiente for es para hacer el calculo de los demás valores desde el final hasta el inicio
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x


# ====================================================
# Método de diferencias finitas para resolver EDOs
# de segundo orden de la forma:
# y'' + p(x)y' + q(x)y = r(x), con condiciones de frontera
# ====================================================

def edo2(p, q, r, h, a, b, y0, yn):
    n = int((b - a) / h) # Número de intervalos 
    x = np.linspace(a, b, n + 1) # Puntos del mallado
    
    # Construimos el sistema lineal de ecuaciones
    A = np.zeros(n - 1)
    B = np.zeros(n - 1)
    C = np.zeros(n - 1)
    D = np.zeros(n - 1)

    for i in range(1, n):
        xi = x[i]
        A[i - 1] = (-h/2)*p(xi)-1  # Subdiagonal A_i 
        B[i - 1] = 2+(h**2)*q(xi)  # Diagonal principal B_i
        C[i - 1] = (h/2)*p(xi)-1  # Superdiagonal C_i
        D[i - 1] = (-h**2)*r(xi)  # Término independiente D_i

    # Ahora vamos a ajustar el vector D para considerar las condiciones de frontera 
    D[0] += ((h/2)*p(x[1])+1) * y0   # Condición en el inicio
    D[n-2] += ((-h/2)*p(x[n])+1) * yn    # Condición en el final 

    # Se resuelve el sistema tridiagonal usando el método de Thomas
    y_inner = metodoDeThomas(A[1:], B, C[:-1], D)

    # Procedemos a hacer la concatenación de las condiciones
    # de frontera al resultado al resultado interno 
    y_total = np.concatenate(([y0], y_inner, [yn]))
    print(y_total)
    return x, y_total

# ==========================
# Definición de las funciones
# ==========================
def p(x): return -1 / x

def q(x): return (1 / (4 * x ** 2)) - 1

def r(x): return 0

# ================================
# Solución exacta para comparar
# ================================
def y_exacta(x):
    return np.sin(6 - x) / (np.sin(5) * np.sqrt(x))


# ===========================
# Prueba del método numérico
# ===========================
h_values = [1, 0.5, 0.2, 0.1, 0.01]

# Graficamos soluciones aproximadas 
plt.figure(figsize=(10, 6))
for h in h_values:
    x, y_aprox = edo2(p, q, r, h, 1, 6, 1, 0)
    plt.plot(x, y_aprox, label=f'h = {h}')

# Solución exacta para comparación
x_exact = np.linspace(1, 6, 500)
plt.plot(x_exact, y_exacta(x_exact), 'k--', label='Solución exacta')

# Configuración del gráfico 
plt.title('Aproximación por Diferencias Finitas vs Solución Exacta')
plt.xlabel('x')
plt.ylabel('y(x)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

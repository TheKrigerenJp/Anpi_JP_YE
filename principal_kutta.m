% Limpiar entorno
clear;
clc;

% Valores para m a probar
m_vals = [10, 20, 50, 100, 250];

% Primera parte: [0, 2] con y(0) = 0
[x1, y1] = runge_kutta_6(0, 2, 0, 100);

% Solución exacta para todo el intervalo [0.01, 10]
x_exact = linspace(0.01, 10, 500);
y_exact = x_exact .* log(x_exact / 2) + 2 * x_exact;

% Gráfica
figure;
hold on;

for m = m_vals
  % Segunda parte: [2, 10] con y(2) = 4
  [x2, y2] = runge_kutta_6(2, 10, 4, m);

  % Unir ambas soluciones
  x_total = [x1, x2];
  y_total = [y1, y2];

  % Graficar cada aproximación
  plot(x_total, y_total, 'DisplayName', ['RK6 m = ', num2str(m)]);
end

% Graficar solución exacta
plot(x_exact, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Exacta');

xlabel('x');
ylabel('y');
title('Soluciones aproximadas vs. solución exacta');
legend('show');
grid on;


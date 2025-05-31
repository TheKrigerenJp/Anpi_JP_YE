function runge_kutta_6()
  clc;
  clear;

  % Parámetros iniciales
  a = 2;
  b = 10;
  y0 = 4;
  ms = [10, 20, 50, 100, 250];

  % Solución exacta
  y_exacta = @(x) log(x / 2) + 2 * x;

  % Crear gráfico
  figure;
  hold on;
  colors = lines(length(ms));

  for i = 1:length(ms)
    m = ms(i);
    [x, y] = runge_kutta_6_2(a, b, y0, m);
    plot(x, y, 'Color', colors(i,:), 'DisplayName', sprintf('RK6 m = %d', m));
  end

  % Graficar la solución exacta
  x_exacto = linspace(a, b, 1000);
  y_exacto = y_exacta(x_exacto);
  plot(x_exacto, y_exacto, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Exacta');

  xlabel('x');
  ylabel('y');
  title('Soluciones aproximadas vs solución exacta');
  legend show;
  grid on;
  hold off;
end

% Función principal Runge-Kutta orden 6
function [x, y] = runge_kutta_6_2(a, b, y0, m)
  h = (b - a) / (m - 1);       % Paso
  x = linspace(a, b, m);      % Vector de puntos x
  y = zeros(1, m);            % Vector de soluciones y
  y(1) = y0;

  for n = 1:(m - 1)
    k1 = h * f(x(n), y(n));
    k2 = h * f(x(n) + h/3, y(n) + k1/3);
    k3 = h * f(x(n) + 2*h/5, y(n) + (4*k1 + 6*k2)/25);
    k4 = h * f(x(n) + h/2, y(n) + (k1 - 12*k2 + 15*k3)/4);
    k5 = h * f(x(n) + 2*h/3, y(n) + (6*k1 + 90*k2 - 50*k3 + 8*k4)/81);
    k6 = h * f(x(n) + 4*h/5, y(n) + (6*k1 + 36*k2 + 10*k3 + 8*k4)/75);
    y(n + 1) = y(n) + (23*k1 + 125*k2 - 81*k3 + 125*k6)/192;
  end
end

% Función f(x, y) = (x + y) / x
function dy = f(x, y)
  dy = (x + y) / x;
end


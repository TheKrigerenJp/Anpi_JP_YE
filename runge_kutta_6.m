function [x, y] = runge_kutta_6(a, b, y0, m)
  % Paso
  h = (b - a) / (m - 1);

  % Vector de puntos en x
  x = linspace(a, b, m);
  y = zeros(1, m);
  y(1) = y0;

  % Iteración con Runge-Kutta de orden 6
  for n = 2:m
    xn = x(n-1);
    yn = y(n-1);

    k1 = h * f(xn, yn);
    k2 = h * f(xn + h/3, yn + k1/3);
    k3 = h * f(xn + 2*h/5, yn + (4*k1 + 6*k2)/25);
    k4 = h * f(xn + h, yn + (k1 - 12*k2 + 15*k3)/4);
    k5 = h * f(xn + 2*h/3, yn + (6*k1 + 90*k2 - 50*k3 + 8*k4)/81);
    k6 = h * f(xn + 4*h/5, yn + (6*k1 + 36*k2 + 10*k3 + 8*k4)/75);

    y(n) = yn + (23*k1 + 125*k2 - 81*k5 + 125*k6)/192;
  end
end

function dy = f(x, y)
  % EDO: y' = x + y
  dy = x + y;
end


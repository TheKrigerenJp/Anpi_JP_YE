function pregunta_1()
  clc;
  clear;

  % ===============================
  % Parámetros del problema
  % ===============================
  a = 2;   % Limite inferior del intervalo
  b = 10;   % Limite superior del intervalo 
  y0 = 4;   % Condicion inicial y(a) = y0
  ms = [10, 20, 50, 100, 250];   % Diferentes valores de particion para probar precision

  % ===============================
  % Solución exacta del problema
  % ===============================
  y_exacta = @(x) log(x / 2) + 2 * x;

  % ===============================
  % Crear la figura para graficar
  % ===============================
  figure;
  hold on;
  colors = lines(length(ms));

  % ===============================
  % Ejecutar Runge-Kutta de orden 6 para cada m con el siguiente for
  % ===============================
  for i = 1:length(ms)
    m = ms(i);
    [x, y] = runge_kutta_6(a, b, y0, m); % Aqui se hace la aproximacion numerica
    plot(x, y, 'Color', colors(i,:), 'DisplayName', sprintf('RK6 m = %d', m));
  end


  % ===============================
  % Graficar la solución exacta
  % ===============================
  x_exacto = linspace(a, b, 1000);
  y_exacto = y_exacta(x_exacto);
  plot(x_exacto, y_exacto, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Exacta');

  % ===============================
  % Ajustes del gráfico
  % ===============================
  xlabel('x');
  ylabel('y');
  title('Soluciones aproximadas vs solución exacta');
  legend show;
  grid on;
  hold off;
end

% ==========================================================
% Método de Runge-Kutta de orden 6 para EDO de primer orden
% ==========================================================
function [x, y] = runge_kutta_6(a, b, y0, m)
  h = (b - a) / (m - 1);      % Paso de integracion 
  x = linspace(a, b, m);      % Vector de puntos x equidistantes 
  y = zeros(1, m);            % Inicializar vector de soluciones 
  y(1) = y0;                  % Asignamos la condicion incial 

  % En el siguiente for se realiza el metodo Runge-Kutta de orden 6
  for n = 1:(m - 1)

    % Calculamos las constantes k segun cada formula
    k1 = h * f(x(n), y(n));
    k2 = h * f(x(n) + h/3, y(n) + k1/3);
    k3 = h * f(x(n) + 2*h/5, y(n) + (4*k1 + 6*k2)/25);
    k4 = h * f(x(n) + h/2, y(n) + (k1 - 12*k2 + 15*k3)/4);
    k5 = h * f(x(n) + 2*h/3, y(n) + (6*k1 + 90*k2 - 50*k3 + 8*k4)/81);
    k6 = h * f(x(n) + 4*h/5, y(n) + (6*k1 + 36*k2 + 10*k3 + 8*k4)/75);

    % Calculamos el valor siguiente de y usando combinacion lineal de las k 
    y(n + 1) = y(n) + (23*k1 + 125*k2 - 81*k3 + 125*k6)/192;
  end
end

% ==========================================================
% Función derivada del problema (definida por el enunciado)
% f(x, y) = (x + y) / x
% ==========================================================
function dy = f(x, y)
  dy = (x + y) / x;
end


clc; clear; close all;
n = 8; % Número de nodos

% Generar una matriz de distancias aleatorias entre 1 y 10
D = [0	3	7	5	2	8	6	9;...
3	0	4	6	1	5	7	10;...
7	4	0	2	9	3	8	5;...
5	6	2	0	7	4	1	3;...
2	1	9	7	0	6	5	8;...
8	5	3	4	6	0	2	7;...
6	7	8	1	5	2	0	4;...
9	10	5	3	8	7	4	0]
%randi([1, 10], n, n);
D = triu(D, 1); % Mantener solo la parte superior
D = D + D'; % Hacer la matriz simétrica
D(eye(n) == 1) = 0; % La diagonal es cero

% Crear un grafo en MATLAB
G = graph(D, 'upper');

% Obtener las distancias de cada arista
L = D; L(L == 0) = inf; % Evitar ceros en la matriz
max_dist = max(G.Edges.Weight); % Máxima distancia
min_dist = min(G.Edges.Weight(G.Edges.Weight > 0)); % Mínima distancia

% Escalar posiciones de los nodos según las distancias
theta = linspace(0, 2*pi, n+1);
theta(end) = [];
x = zeros(1, n);
y = zeros(1, n);

for i = 2:n
    scale_factor = (G.Edges.Weight(i-1) - min_dist) / (max_dist - min_dist) + 0.5;
    x(i) = cos(theta(i)) * scale_factor;
    y(i) = sin(theta(i)) * scale_factor;
end

% Dibujar el grafo con distancias escaladas
figure;
H = plot(G, 'XData', x, 'YData', y, 'LineWidth', 1.5, 'EdgeLabel', G.Edges.Weight);
%title('Grafo con distancias ajustadas visualmente');
set(gca, 'XColor', 'none', 'YColor', 'none'); % Ocultar ejes

% Definir la matriz de distancia
D = [0	3	7	5	2	8	6	9;
     3	0	4	6	1	5	7	10;
     7	4	0	2	9	3	8	5;
     5	6	2	0	7	4	1	3;
     2	1	9	7	0	6	5	8;
     8	5	3	4	6	0	2	7;
     6	7	8	1	5	2	0	4;
     9	10	5	3	8	7	4	0];

% Convertir la matriz de distancia en coordenadas 2D usando MDS
Y = cmdscale(D);

% Crear un grafo desde la matriz de distancia (asume conexiones completas)
G = graph(D, 'upper', 'omitselfloops');

% Graficar el grafo
figure;
h = plot(G, 'XData', Y(:,1), 'YData', Y(:,2), 'EdgeLabel', G.Edges.Weight);

% Personalizar el gráfico
h.NodeLabel = {'1','2','3','4','5','6','7','8'}; % Etiquetas de nodos
title('Visualización del Grafo con Distancias entre Nodos');
axis equal;
grid on;

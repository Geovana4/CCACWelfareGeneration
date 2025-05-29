% Extraer la matriz de potencias PjiSol en el tiempo
PjiSol_tiempo = y2(:, 1:(generadores * consumidores));
PjiSol_filtro = y(:, 1:(generadores * consumidores));
PjiSol_filtro2 = y3(:, 1:(generadores * consumidores));

% Inicializar vector de bienestar en el tiempo
num_tiempo = length(x2);
num_filtro = length(x);
num_filtro2 = length(x3);
Welfare_time = zeros(num_tiempo, 1);
Welfare_filtro = zeros(num_filtro, 1);
Welfare_filtro2 = zeros(num_filtro2, 1);

% Calcular bienestar en cada instante de tiempo
for t = 1:num_tiempo
    Welfare_time(t) = Welfarejgen(PjiSol_tiempo(t, :), Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator);
end

for t = 1:num_filtro
    Welfare_filtro(t) = Welfarejgen(PjiSol_filtro(t, :), Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator);
end


for t = 1:num_filtro2
    Welfare_filtro2(t) = Welfarejgen(PjiSol_filtro2(t, :), Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator);
end

% Calcular el bienestar de referencia con Pij_opt_new
Welfare_opt = Welfarejgen(Pij_opt_new, Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator);
Welfare_in = Welfarejgen(Xo, Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator);
WelFil = Welfare_filtro(end);
WelFil2 = Welfare_filtro2(end);
WelTime = Welfare_time(end);
disp(Welfare_opt)
disp('PI')
disp(Welfare_in)
disp('Inicial')
disp(WelFil)
disp('Con filtro 0.01')
disp(WelFil2)
disp('Con filtro2 0.001')
disp(WelTime)
disp('Sin filtro')
num_pares = generadores * consumidores;

% Crear la figura
figure(1);
subplot(3,1,1);
hold on;
for i = 1:num_pares
    plot(1:num_tiempo, PjiSol_tiempo(:, i), 'LineWidth', 1);
end
xlabel('Step');
ylabel('$P_{ji}$ (kWh)','Interpreter', 'latex');
title('A');
xlim([0, 200])
grid on; 
hold off;
subplot(3,1,2);
hold on;
for i = 1:num_pares
    plot(1:num_filtro, PjiSol_filtro(:, i), 'LineWidth', 1);
end
xlabel('Step');
ylabel('$P_{ji}$ (kWh)','Interpreter', 'latex');
title('B');
xlim([0 200])
legend({'$\tau=0.01$'}, 'Interpreter', 'latex',...
       'Location', 'east'); % Corregido símbolo tau en LaTeX
grid on;
hold off;
subplot(3,1,3);
hold on;
for i = 1:num_pares
    plot(1:num_filtro2, PjiSol_filtro2(:, i), 'LineWidth', 1);
end
xlabel('Step');
ylabel('$P_{ji}$ (kWh)','Interpreter', 'latex');
title('C');
xlim([0 200])
legend({'$\tau=0.001$'}, 'Interpreter', 'latex',...
       'Location', 'east'); % Corregido símbolo tau en LaTeX
grid on;
f = gcf;
exportgraphics(f,'pot144.png','Resolution',300)
hold off;


% Reestructurar PjiSol_tiempo y Pij_opt_new en matrices de tamaño (generadores x consumidores)
PjiSol_tiempoM = reshape(PjiSol_filtro2, [num_filtro2, generadores, consumidores]);
Pij_opt_newM = reshape(Pij_opt_new, [generadores, consumidores]);

% Cálculo del consumo total en el tiempo (suma por filas, es decir, sobre generadores)
consumoF_tiempo = squeeze(sum(PjiSol_tiempoM, 2)); % Suma sobre la segunda dimensión (generadores)
consumoF_ref = sum(Pij_opt_newM, 1);  % Suma sobre filas (generadores) para referencia

% Cálculo de la generación total en el tiempo (suma por columnas, es decir, sobre consumidores)
generacionF_tiempo = squeeze(sum(PjiSol_tiempoM, 3)); % Suma sobre la tercera dimensión (consumidores)
generacionF_ref = sum(Pij_opt_newM, 2);  % Suma sobre columnas (consumidores) para referencia

figure(2);
plot(1:num_tiempo, Welfare_time, 'b', 'LineWidth', 1.5); hold on;
plot(1:num_filtro, Welfare_filtro, 'g', 'LineWidth', 1.5); hold on;
plot(1:num_filtro2, Welfare_filtro2, 'r', 'LineWidth', 1.5); hold on;
yline(Welfare_opt, 'k--', 'LineWidth', 1.3); % Línea de referencia
xlabel('Step');
xlim([0 200])
ylabel('$W_j$', 'Interpreter', 'latex', 'Color', 'k'); 
f = gcf;
%title('A');
legend('RD','RD $\tau=0.01$', 'RD $\tau=0.001$','PI','Interpreter', 'latex');
exportgraphics(f,'welfare144.png','Resolution',300)
grid on;

figure(3);
x0=10;
y0=10;
width=900;
height=700
set(gcf,'position',[x0,y0,width,height])

% Configurar colores únicos para consumidores y generadores
colores_consumo = parula(consumidores); % Paleta de colores para consumidores
colores_generacion = parula(generadores); % Paleta de colores para generadores

% Gráfica 1: Consumo en el tiempo
subplot(2,1,1);
hold on;

% Graficar cada consumidor con su color y almacenar handles
h_consumo = zeros(1, consumidores);
h_ref = zeros(1, consumidores);
Di_ref = repmat(Di, num_filtro2, 1);

for i = 1:consumidores
    h_consumo(i) = plot(1:num_filtro2, consumoF_tiempo(:,i),...
                       'Color', colores_consumo(i,:),...
                       'LineWidth', 1.2,Marker='|');
    h_ref(i) = plot(1:num_filtro2, Di_ref(:,i),...
                   'Color', colores_consumo(i,:),...
                   'LineWidth', 1.2, 'LineStyle', '--'); % Línea discontinua para diferenciar ref
end

% Línea de referencia (PI)
h_ref_consumo = yline(consumoF_ref, 'k:', 'LineWidth', 1.3);
xlim([0 40])
% Configurar etiquetas y leyenda
xlabel('Step');
ylabel('$\sum_{i \in \mathcal{I}} P_{ji}$ (kWh)', 'Interpreter', 'latex');
title('A');
%xlim([0 num_filtro2])

% Crear etiquetas dinámicas: D1, ref1, D2, ref2,..., PI
leyenda_consumo = [reshape([arrayfun(@(x) sprintf('D_%d', x), 1:consumidores, 'UniformOutput', false); ...
                           arrayfun(@(x) sprintf('ref_%d', x), 1:consumidores, 'UniformOutput', false)], 1, []), {'PI'}];

legend(leyenda_consumo, 'Location', 'east');

grid on;
hold off;

% Gráfica 2: Generación en el tiempo
subplot(2,1,2);
hold on;

% Graficar cada generador con su color y almacenar handles
h_generacion = zeros(1, generadores);
h_ref2 = zeros(1, generadores);
Gj_ref = repmat(Gj, num_filtro2, 1);

for i = 1:generadores
    h_generacion(i) = plot(1:num_filtro2, generacionF_tiempo(:,i),...
                         'Color', colores_generacion(i,:),...
                         'LineWidth', 1.2,Marker='|');
    h_ref2(i) = plot(1:num_filtro2, Gj_ref(:,i),...
                   'Color', colores_generacion(i,:),...
                   'LineWidth', 1.2, 'LineStyle', '--'); % Línea discontinua para diferenciar ref
end

% Línea de referencia (PI)
h_ref_generacion = yline(generacionF_ref, 'k:', 'LineWidth', 1.3);

% Configurar etiquetas y leyenda
xlabel('Step');
ylabel('$\sum_{j \in \mathcal{J}} P_{ji}$ (kWh)', 'Interpreter', 'latex');
title('B');
% Crear etiquetas dinámicas: G1, ref1, G2, ref2,..., PI
leyenda_generacion = [reshape([arrayfun(@(x) sprintf('G_%d', x), 1:generadores, 'UniformOutput', false); ...
                              arrayfun(@(x) sprintf('ref_%d', x), 1:generadores, 'UniformOutput', false)], 1, []), {'PI'}];

legend(leyenda_generacion, 'Location', 'east');
xlim([0 40]);
grid on;
f = gcf;
exportgraphics(f, 'restri144.png', 'Resolution', 300);
hold off;


function Wj_total = Welfarejgen(y, Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator)
    Pij = reshape(y, [generadores, consumidores]);
    %disp(Pij)
    Pgb = 114;
%     Wj = arrayfun(@(j) lamdaj(j) * Dj(j) - thetaj(j) * Dj(j)^2 + ...
%         Pgb.* log(sum(Pij(j, :) ./(pii+1))) - ...
%         a(j) * (sum(Pij(j, :))^2) - ...
%         b(j) * sum(Pij(j, :)), generator);
     Wj = arrayfun(@(j) lamdaj(j) * Dj(j) - thetaj(j) * Dj(j)^2 + ...
         sum(Pij(j, :) .*pii) - ...
         a(j) * (sum(Pij(j, :))^2) - ...
         b(j) * sum(Pij(j, :)), generator);
    %disp(Wj)
%     Reco = log(Pij./(pii+1))
%     Wj = arrayfun(@(j) lamdaj(j) * Dj(j) - thetaj(j) * Dj(j)^2 + ...
%         Pgb.*sum(Reco(j,:)) - ... 
%         a(j) * (sum(Pij(j, :))^2) - ...
%         b(j) * sum(Pij(j, :)), generator);
    Wj_total = -sum(Wj);
end
%%%Nuevo codigo con filtro pata generadores 
clc
%clear all
% Cargar el archivo Excel
file_path = 'Demandaóptima_Comunidad8pCon.xlsx';
file_path2 = 'Precios.xlsx';
df = readtable(file_path);
df2 = readtable(file_path2);
% Extraer las columnas de interés
G1 = df.Glim1;
G2 = df.Glim2;
G3 = df.Glim3;
G4 = df.Glim4;
G5 = df.Glim5;
G6 = df.Glim6;
G7 = df.Glim7;
G8 = df.Glim8
D1 = df.D11; 
D2 = df.D22;
D3 = df.D33;
D4 = df.D44;
D5 = df.D55;
D6 = df.D66;
D7 = df.D77;
D8 = df.D88;

pii1 = df2.pii1;
pii2 = df2.pii2;
pii3 = df2.pii3;
pii4 = df2.pii4;
pii5 = df2.pii5;
pii6 = df2.pii6;
pii7 = df2.pii7;
pii8 = df2.pii8;
piiM = [pii1,pii2,pii3,pii4,pii5,pii6,pii7, pii8]
% Crear matrices de generación y demanda
generation = [G1, G2, G3, G4, G5, G6, G7,G8];
demand = [D1, D2, D3, D4, D5, D6,D7,D8];

% Parámetros de entrada
theta = [1/2, 1/2, 1/2, 1/2, 1/2, 1/2,1/2,1/2];
lamda = [100, 100, 100, 100, 100, 100,100,100];
etha0 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1,0.1];
Pgs = 1250;
Pgb = 114;
precioin = [Pgb, Pgb, Pgb, Pgb, Pgb, Pgb,Pgb, Pgb ];
pii0_G = [Pgb, Pgb, Pgb, Pgb, Pgb, Pgb, Pgb, Pgb];

% Número de prosumers
prosumers = 8;
N = 1:prosumers;
PmasCP = 1:prosumers;

% Coeficientes a y b
a= [4*0.089,0.069,0.067,0,0,0,0,0]; %Termica, Sol+W, Sol+W, eol, eol 10 veces los costos iniciales 
a= [6.0865*a];

b= [3.93*52,32,28,47,37,0,0,0];
b = [6.0865*b];

% Seleccionar un tiempo específico (t=6)
t =14;
Glim = generation(t, :);
Dopt = demand(t, :);
piif = piiM(t,:)
% Filtrar datos según la condición Glim[n] / Dopt[n] <= 1
Di = Dopt(Glim ./ Dopt <= 1) - Glim(Glim ./ Dopt <= 1);
Gi = Glim(Glim ./ Dopt <= 1);
thetai = theta(Glim ./ Dopt <= 1);
lamdai = lamda(Glim ./ Dopt <= 1);
etha = etha0(Glim ./ Dopt <= 1);

% Filtrar datos según la condición Glim[n] / Dopt[n] >= 1
Dj = Dopt(Glim ./ Dopt >= 1);
Gj = Glim(Glim ./ Dopt >= 1) - Dopt(Glim ./ Dopt >= 1);
GjBound2 = Glim(Glim ./ Dopt >= 1);
thetaj = theta(Glim ./ Dopt >= 1);
lamdaj = lamda(Glim ./ Dopt >= 1);

% Número de consumidores y generadores
consumidores = length(Di);
consumer = 1:consumidores;
generadores = length(Gj);
generator = 1:generadores;

% Filtrar coeficientes a y b para generadores
a = a(Glim ./ Dopt >= 1);
b = b(Glim ./ Dopt >= 1);
pii = piif(Glim ./ Dopt < 1)

%%%%%%%%%%%%%%%%%%Condiciones iniciales
if sum(Gj) > sum(Di)
    simplex = sum(Di);
    % Crear la matriz Xo
    Xo = arrayfun(@(i) Di(i) / generadores, 1:consumidores, 'UniformOutput', false);
    Xo = repmat(Xo, generadores, 1);
    Xo = cell2mat(Xo);
    Xo = Xo(:)'; % Aplanar la matriz
    Xo2 = 0.1 * ones(1, generadores + consumidores);
    Xoo = [Xo, Xo2];
    disp(Xoo);
else
    simplex = sum(Gj);
    % Crear la matriz Xo
    Xo = arrayfun(@(j) Gj(j) / consumidores, 1:generadores, 'UniformOutput', false);
    Xo = repmat(Xo, consumidores, 1);
    Xo = cell2mat(Xo);
    Xo = Xo'; % Transponer y aplanar
    Xo = Xo(:)';
    Xo2 = ones(1, generadores + consumidores);
    Xoo = [Xo, Xo2];
    disp(Xoo);
end

PijSol = reshape(Xo, [generadores, consumidores]);
%pii= [268.71,292.65,280.25,344.74] %14pm
%pii= [335.37,258.36, 291.23, 281.15, 348.51]
t_span = [0,1]; % Intervalo de tiempo
t_eval = linspace(t_span(1), t_span(2), 1000); % Puntos de evaluación
%t_eval = [0:0.001:1]
y0 = Xoo; % Condición inicial
disp(y0); % Imprimir valores iniciales

timein = tic; % Iniciar temporizador

% Resolver la ecuación diferencial
options = odeset('RelTol',1e-10, 'AbsTol',1e-10);
x0_filtros = zeros(generadores + consumidores, 1); % Un filtro por generador y consumidor

% Estado inicial combinado
x0_combined = y0;

%sol2 = ode15s(@(t, X) ReplicadorWjSol2(t, X, a, b, solutionp, simplex, generadores, consumidores, generator, consumer, Di, Gj), ...
            %  t_span, y0, options);

tau =0.1; % Constante de tiempo del filtro
% Llamar al solver
[x2,y2] = ode23s(@(t, y) ReplicadorWjSol2(t, y, a, b, pii, simplex, generadores, consumidores, generator, consumer, Di, Gj, tau), t_eval, x0_combined);

timefin_sinFil = toc(timein); % Finalizar temporizador
disp(timefin_sinFil); % Mostrar el tiempo de ejecución


function totales = ReplicadorWjSol2(t, y, a, b, pii, simplex, generadores, consumidores, generator, consumer, Di, Gj, tau)
    % Extraer parámetros desde la estructura
    VelGrad = 1e5;
    bgrande = 1e5;
    VelRD = 0.1;
    Pgs = 1250; 

    % Extraer variables de estado
    total_Pij = generadores * consumidores;
    Pij = y(1:total_Pij);
    Pij = reshape(Pij, [generadores, consumidores]);

    lamdaub = y(total_Pij + 1 : total_Pij + generadores);
    betaub = y(total_Pij + generadores + 1 : total_Pij + generadores + consumidores);

    %lamdaub_filt = y(total_Pij + generadores + consumidores + 1 : total_Pij + 2*generadores + consumidores);
    %betaub_filt = y(total_Pij + 2*generadores + consumidores + 1 : total_Pij + 2*generadores + 2*consumidores);
    H = arrayfun(@(j) 2 * a(j) * sum(Pij(j, :)) + b(j), generator);
    %R = arrayfun(@(i) pii(i), consumer);
    Pgb = 114;
    %R = Pgb.*((1 ./ Pij) * ones(consumidores, 1));
    R = arrayfun(@(i) pii(i), consumer);
    % Calcular dPji
    dPji = arrayfun(@(j, i) pii(i) - H(j)- lamdaub(j) - betaub(i) + bgrande, ...
        repmat((1:generadores)', 1, consumidores), repmat(1:consumidores, generadores, 1));

    % Calcular dPjimean
    dPjimean = sum(sum(Pij .* dPji)) / simplex;

    % Calcular dPijP
    dPij = Pij .* VelRD .* (dPji - dPjimean);
    %dPijP = dPijP(:)

    % Calcular Gbetaub y Glamdaub (asegurar que sean vectores columna)
    Gbetaub = arrayfun(@(i) VelGrad * betaub(i) * (sum(Pij(:, i)) - Di(i)) + 100, consumer);
    Gbetaub = Gbetaub(:); % Convertir a columna

    Glamdaub = arrayfun(@(j) VelGrad * lamdaub(j) * (sum(Pij(j, :)) - Gj(j)) +100, generator);
    Glamdaub = Glamdaub(:); % Convertir a columna


    % Derivadas de los filtros para lambda y beta
%     dGlamdaub_filt = (Glamdaub-lamdaub_filt) / tau;
%     dGbetaub_filt = (Gbetaub-betaub_filt) / tau;
    totales = [dPij(:); Glamdaub; Gbetaub];
end
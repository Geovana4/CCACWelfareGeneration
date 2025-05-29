
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
lb_Pij = zeros(numel(Xo), 1);
ub_Pij = inf(numel(Xo), 1);

%pii= [268.71,292.65,280.25,344.74]
%pii= [335.37,258.36, 291.23, 281.15, 348.51]
y0 = Xoo; % Condición inicial
disp(y0); % Imprimir valores iniciales

timein = tic; % Iniciar temporizador
[Pij_opt_new, welfarej_optj, tiempo1] = WjOpt(Xo, lb_Pij, ub_Pij, Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator,Gj,Di);
timefin = toc(timein); % Finalizar temporizador
disp(timefin); % Mostrar el tiempo de ejecución
function [c, ceq] = ConstraintsWj(x, generadores, consumidores, Gj,Di)
    c = [];
    ceq = [];
    Pij = reshape(x, [generadores, consumidores]);
    % Restricciones de capacidad
    c = [c; sum(Pij,2) - Gj']; % sum(Pij) <= Gj
    % Restricciones de demanda
    c = [c; sum(Pij,1)'-Di']; % sum(Pij)<= Di
    % Restricciones de igualdad
    if sum(Di) <= sum(Gj)
        ceq = sum(Pij,1)' - Di';
    else
        ceq = sum(Pij,2) - Gj';
    end
end

function [Wj_opt, welfarej_optj, tiempo] = WjOpt(sli0, lb, ub, Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator,Gj,Di)
    tic;
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', 'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-10);
    [Wj_opt, fval] = fmincon(@(x)Welfarejgen(x, Dj, a, b, thetaj, lamdaj, pii, generadores, consumidores, generator), sli0, [], [], [], [], lb, ub, @(x)ConstraintsWj(x, generadores, consumidores, Gj,Di), options);
    %Pij_opt = Wj_opt;
    welfarej_optj = -fval;
    tiempo = toc;
end

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

function [Costos, re, pagos] = Costosk(Pij, a, b, pii, generator)
    Costos = arrayfun(@(j) a(j)*sum(Pij(j,:))^2 + b(j)*sum(Pij(j,:)), generator);
    re = sum(Pij .* pii, 2)';
    pagos = pii.*sum(Pij,1);
end

function [demandafinal, Di, generacionfinal, Gj] = PruebaRestriG(Di, Pij, Gj)
    demandafinal = sum(Pij,1);
    generacionfinal = sum(Pij,2)';
end
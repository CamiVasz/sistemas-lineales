
%===== Archivo de ejemplo para la linealización y obtención de funciones de
%transferencia utilizando MATLAB=====
%NOTA: Se debe tener instalado el Toolbox Symbolic Math



%====Primero se declaran las variables del modelo como simbolicas====
syms x c ce u

%==== Se definen los parámetros del modelo 

syms a n b h r0 r1 p k0 M mu la
% a  = 0.01;
% n  = 0.01;
% b  = 0.2;
% h  = 0.2;
% r0 = 0.1;
% r1 = 0.01;
% p  = 2;
% k0 = 5;
% M  = 0.5;
% mu = 0.1;
% la = 0.01;


%==== Se escriben las ecuaciones del modelo ====
dx  = x*((r0 + r1*(c/x))/((c/x)+1))*(1-(x/((p*(c/x) + k0)/((c/x) + 1))));
dce = -a*ce*x + b*c - h*ce + u;
dc  = n*ce*x - mu*c;

%==== Se crea el vector de ecuaciones ====
FF = [dx dce dc]; % Vector de Ecuaciones 

%==== Se define el vector de salidas ====

Y= x; %Ecuación de la salida

%==== Se declaran las variables de estado y entradas ====
x = [x ce c]; % Variables de estado- Variables que se encuentran en las derivadas del modelo
U = [u];  % Vector de
entradas

%==== Se realiza la derivadas parciales por medio del comando
%jacobian() del vector de ecuaciones con respecto de las variables de estado y las entradas====
A = jacobian(FF,x);
BM = jacobian(FF,U);
C = jacobian(Y,x);
D = jacobian(Y,u);

%==== Se define el punto de operación al rededor del cual se realizará
%lalinealización ====
u = 4;
x = 2.9;
c = 6.8;
ce = 23.4;

%==== Se realiza la lineación en el punto de operación ====
A = eval(A);
BM = eval(BM);
C = eval(C);
D = eval(D);

[As, BMs, Cs, Ds] = linmodv5('biomasaLineal', [x ce c], 4)

plot(tsim, xnl, tsim, xli)
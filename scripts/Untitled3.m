
%===== Archivo de ejemplo para la linealizaci�n y obtenci�n de funciones de
%transferencia utilizando MATLAB=====
%NOTA: Se debe tener instalado el Toolbox Symbolic Math



%====Primero se declaran las variables del modelo como simbolicas====
syms x c ce u

%==== Se definen los par�metros del modelo 

syms a n b h r0 r1 p k0 M mu la

n  = 0.01;
b  = 0.2;
h  = 0.2;
r0 = 0.1;
r1 = 0.01;
p  = 2;
k0 = 5;
M  = 0.5;
mu = 0.1;
la = 0.01;


%==== Se escriben las ecuaciones del modelo ====
dx  = x*((r0 + r1*(c/x))/((c/x)+1))*(1-(x/((p*(c/x) + k0)/((c/x) + 1))));
dce = -a*ce*x + b*c - h*ce + u;
dc  = n*ce*x - mu*c;

%==== Se crea el vector de ecuaciones ====
FF = [dx dce dc]; % Vector de Ecuaciones 

%==== Se define el vector de salidas ====

Y= x; %Ecuaci�n de la salida

%==== Se declaran las variables de estado y entradas ====
x = [x ce c]; % Variables de estado- Variables que se encuentran en las derivadas del modelo
U = [u];  % Vector de entradas

%==== Se realiza la derivadas parciales por medio del comando
%jacobian() del vector de ecuaciones con respecto de las variables de estado y las entradas====
A = jacobian(FF,x);
BM = jacobian(FF,U);
C = jacobian(Y,x);
D = jacobian(Y,u);

%==== Se define el punto de operaci�n al rededor del cual se realizar�
%lalinealizaci�n ====
u = 4;
x = 2.9;
c = 6.8;
ce = 23.4;

%==== Se realiza la lineaci�n en el punto de operaci�n ====
A = eval(A);
BM = eval(BM);
C = eval(C);
D = eval(D);

a  = 0.01;
A
sys = ss(A,BM,C,D)
G = tf(sys)
GD = c2d(G,0.5)

figure
step(sys)
hold on
step(G)
hold on 
step(GD)

figure
impulse(sys)
hold on
impulse(G)
hold on 
impulse(GD)

T = 0:0.5:180;
uu = sin(T) + 4;
figure
lsim(sys, uu, T)
hold on
lsim(G, uu, T)
hold on 
lsim(GD, uu, T)

% Funcion de ponderacion discreta
syms z
zz = [z z^2 z^3 z^4];
[num, den] = tfdata(GD, 'v');
iz = (zz * num')/(zz * den');
% transformada inversa z
ff = iztrans(iz);
fpd = zeros(10,1);
% Se hallan los primeros 10 términos
for i = 1:11
    n = sym(i-1);
    n0 = eval(ff);
    fdp(i) = eval(n0);
end

fdp = [-0.000001416533741415915,
    -0.00000825786463191854,
    -0.00002384029688964825,
    -0.000048047247438898,
    -0.00007930902293923858,
    -0.0001162670126816055,
    -0.0001577471214501992,
    -0.0002027364960132934,
    -0.0002503631382311275,
    -0.0002998780480342695,
    -0.0003506395835840654];
T = 0:0.5:5;
plot(0:0.25:2.5, fdp, 'o')
hold on
%figure
[x,y] = impulse(GD,T)


% reduccion
%syms s
[num, den] = tfdata(G, 'v');
G0 = num(end)/den(end);
den = (s-0.0474 - 0.0156i)*(s-0.0474 + 0.0156i);
dr = flip(eval(coeffs(den)));
%K = G0 / (num(end) / den);
K = (-9.89*10^(-5)*dr(end))/0.000655;
aa = pole(G)
aa = aa(2:end);
rs = tf(zpk([], aa, K));
step(rs);
hold on
step(G)
hold off

T = 0:0.5:180;
uu = sin(T) + 4;
figure
lsim(rs, uu, T)
hold on
lsim(G, uu, T)


GR = tf(reduce(G,2));
step(GR);
hold on
step(G)
hold off

T = 0:0.5:180;
uu = sin(T) + 4;
figure
lsim(GR, uu, T)
hold on
lsim(G, uu, T)

Gz = (-1.971*10^(-6)*z^2 - 7.54*10^(-6)*z - 1.803*10^(-6))/(z^3 - 2.83*z^2 + 2.666*z - 0.8361)
iz = iztrans(Gz)

fpd4 = zeros(10,1);
for i = 1:11
    n = sym(i-1);
    n0 = eval(iz);
    fdp4(i) = eval(n0);
end
fdp2 = [-3.7947e-19,
    -1.9710e-06,
    -1.3118e-05,
    -3.3672e-05,
    -6.1967e-05,
    -9.6566e-05,
    -1.3623e-04,
    -1.7990e-04,
    -2.2666e-04,
    -2.7574e-04,
    -3.2648e-04]
fdp3 = [0.,-1.971*10^-6,-0.0000131179,-0.0000336721,-0.0000619675,-0.0000965661,-0.00013623,-0.000179897,-0.000226658,-0.000275738,-0.000326481];
T = 0:0.5:5;
plot(T, fdp4, 'o')
hold on
%figure
impulse(GD,T)

% Lugar de las raices
s = tf('s');
GRL = (2.9*s^2 +0.4*s + 0.0008)/(s^3 + 0.3*s^2 + 0.028*s + 0.0006)
rlocus(GRL)

%Diagrama de bode
bode(G)
% modelo reducido
bode(rs)

%modelo discretizado
bode(GD)

T = 0:0.5:180;
uu = sin(0.1*T);
lsim(GD, uu, T)

margin(GD)
% periodo de muestreo
% ancho de banda
wb = 0.034;
Ts2 = [pi/(15*wb), pi/(2.5*wb)];

step(sys)
vf = -0.151;
dd = 0.1*vf;
nn = 0.9*vf;

Tr = 76.8 - 14.6;

Ts1 = [Tr/10, Tr/2];

Ts = 2*wb;

% control PID

% Curva de reaccion
x1 = 25.27;
x2 = 36.12;
y1 = -0.04505;
y2 = -0.07423;

m = (y2 - y1)/(x2 - x1);

b = y2 - m*x2;

L = -b/m;

% ======================================================
R = m;
L = L + 0.5 / 2;

% Metodo Ziegler Nichols

ZN = [  1/(R*L),   Inf,     0;
      0.9/(R*L), 3*L,     0;
      1.2/(R*L), 2*L, 0.5*L];
  
% Metodo de Chien - Hrones - Reswick

CHR = [ 0.3/(R*L),   Inf,      0, 0.7/(R*L),     Inf,      0;
        0.6/(R*L),   4*L,      0, 0.7/(R*L),   2.3*L,      0;
       0.95/(R*L), 2.4*L, 0.42*L, 1.2/(R*L),     2*L, 0.42*L];
 
% periodo de muestreo
T = 0.5;

% P
% Parámetros ZN
Kp = ZN(1,1);
Ti = ZN(1,2);
Td = ZN(1,3);

% CHR
Kp = CHR(1,1);
Ti = CHR(1,2);
Td = CHR(1,3);

% PI
% Parámetros ZN
Kp = ZN(2,1);
Ti = ZN(2,2);
Td = ZN(2,3);

% CHR
Kp = CHR(2,1);
Ti = CHR(2,2);
Td = CHR(2,3);

% PID
% parametroz ZN
Kp = ZN(3,1);
Ti = ZN(3,2);
Td = ZN(3,3);

% CHR
Kp = CHR(3,1);
Ti = CHR(3,2);
Td = CHR(3,3);

% coeficientes funcion de transferencia
q0 = Kp*(1+(T/(2*Ti)) + Td/T);
q1 = -Kp*(1-(T/(2*Ti))+((2*Td)/T));
q2 = (Kp*Td)/T;

% ASIGNACIÓN DE POLOS
% FUNCION DE TRANSFERENCIA DE LA PLANTA

syms q0 q1 q2 r1 r2 z
%r2 = -r1-1;
Bp = (-1.971e-06*z^2 - 7.542e-06*z - 1.803e-06);
Ap = (z^3 - 2.83*z^2 + 2.666*z - 0.8361);
Qc = q0*z + q1;
Pc = (z-r1)*(z-1);

s1 = factor(Ap,'FactorMode','real');
s2 = factor(Bp,'FactorMode','real');

Ad = (z - 0.4)*(z - 0.5)*(z - 0.3);
pol = s1(1)*Pc + s2(1)*s2(3)*Qc - Ad;
C = coeffs(pol,z);
S = solve(C==[0 0 0]);

% Analisis de controlabilidad

Mc = ctrb(A,BM);
rank(Mc)
cond(Mc)

% Realimentacion de estado

% sistema discreto
sysd = c2d(sys, 0.5);
A1 = sysd.A;
B1 = sysd.B;
poles = [0,0.4,0.1];
K = place(A1, B1, [0,0.4,0.1]);

% control del error en estado estacionario
poles1 = [0,0.1,0.2,0.3];
[K, L] = poles_int(sysd, poles1)
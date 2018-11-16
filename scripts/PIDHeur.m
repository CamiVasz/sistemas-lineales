% Curva de reaccion
x1 = 1.4;
x2 = 1.31;
y1 = 1.01;
y2 = 0.966;

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
   
   
   
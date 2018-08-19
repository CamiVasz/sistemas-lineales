% METODO DE EULER
function [ti, xi, cei, ci] = eu(s)

syms x ce c cb

% parametros

a  = 0.01;
n  = 0.01;
b  = 0.2;
h  = 0.2;
r0 = 0.1;
r1 = 0.01;
p  = 2;
k0 = 5;
mu = 0.1;

% definimos el tamaño del paso
t0 = 0;
F  = 200;

%entrada
u   = 0.12;

% ecuaciones diferenciales
dx  =  x*((r0 + r1*cb)/(cb+1))*(1-(x/((p*cb + k0)/(cb + 1))));
dce = -a*ce*x + b*c - h*ce + u;
dc  =  n*ce*x - mu*c;

% Condiciones iniciales
x0  = 2;
c0  = 0;
ce0 = 2;

% tiempo de simulacion
ti      = t0 : s : F;
% vectores de solucion
xi      = zeros(1, length(ti));
ci      = zeros(1, length(ti));
cei     = zeros(1, length(ti));
xi(1)   = x0;
ci(1)   = c0;
cei(1)  = ce0;

for i = 2 : length(ti)
    t     = ti(i-1);
    x     = xi(i-1);
    ce    = cei(i-1);
    c     = ci(i-1);
    if xi(i-1) == 0
        cb = 0;
    else
        cb = ci(i-1)/xi(i-1);
    end
    xi(i) = xi(i-1) + s*eval(dx);
    ci(i) = ci(i-1) + s*eval(dc);
    cei(i) = cei(i-1) + s*eval(dce);
end
end
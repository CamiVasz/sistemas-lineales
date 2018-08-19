%Runge-Kutta
function [tr, xr, cer, cr] = rungeKutta(ts)
syms x ce c cb
% definimos el tamaño del paso
t0 = 0;
F  = 200;

% tiempo de simulacion
tr  = t0 : ts : F;
n   = length(tr);
xr  = zeros(1,n);
cer = zeros(1,n);
cr  = zeros(1,n);

% Condiciones iniciales
xr(1)  =  2;
cer(1) =  2;
cr(1)  = 0;

u = 0.12;

% parametros

a  = 0.01;
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


% Ecuaciones
dx  = x*((r0 + r1*cb)/(cb+1))*(1-(x/((p*cb + k0)/(cb + 1))));
dce = -a*ce*x + b*c - h*ce + u;
dc  = n*ce*x - mu*c;


for i = 1 : length(tr)-1
    x     = xr(i);
    ce   = cer(i);
    c     = cr(i);
   if xr(i) == 0
        cb = 0;
    else
        cb = cr(i)/xr(i);
   end
    
   k1 = ts*eval(dx);
   m1 = ts*eval(dce);
   n1 = ts*eval(dc);
   
   x     = xr(i) + k1/2 ;
   ce   = cer(i) + m1/2;
   c     = cr(i) + n1/2;
   
   k2 = ts*eval(dx);
   m2 = ts*eval(dce);
   n2 = ts*eval(dc); 
   
   x     = xr(i) + k2/2 ;
   ce   = cer(i) + m2/2;
   c     = cr(i) + n2/2;
   
   k3 = ts*eval(dx);
   m3 = ts*eval(dce);
   n3 = ts*eval(dc); 
   
   x     = xr(i) + k3/2 ;
   ce   = cer(i) + m3/2;
   c     = cr(i) + n3/2;
   
   k4 = ts*eval(dx);
   m4 = ts*eval(dce);
   n4 = ts*eval(dc);
   
   xr(i+1) = xr(i) + (k1+(k2*2)+(k3*2)+k4)/6;
   cer(i+1) = cer(i) + (m1+(m2*2)+(m3*2)+m4)/6;
   cr(i+1) = cr(i) + (n1+(n2*2)+(n3*2)+n4)/6;
end
end



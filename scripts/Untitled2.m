% Curva de linealidad

U = 0:0.01:8;
x_cl = zeros(length(U), 1);
c_cl = zeros(length(U), 1);
ce_cl = zeros(length(U), 1);

for i = 1:length(U)
    u = U(i);
    sim('biomasa');
    x_cl(i) = xsim(end);
    c_cl(i) = csim(end);
    ce_cl(i) = cesim(end,1);
end

figure
%subplot(1,3,1)
hold on
%plot(U, U)
plot(U, x_cl)
hold off
subplot(1,3,2)
plot(U, c_cl)
subplot(1,3,3)
plot(U, ce_cl)

% punto de equilibrio

% variables de estado
syms x c ce
% parametros

% entrada a linealizar
u = 4;
n = 0.01;
a = 0.01;
b = 0.2;
h = 0.2;
mu = 0.1; 
r0 = 0.1;
k0 = 5;
r1 = 0.01;
p = 2;

% Ecuaciones
dx  = x*((r0 + r1*(c/x))/((c/x)+1))*(1-(x/((p*(c/x) + k0)/((c/x) + 1))));
dce = -a*ce*x + b*c - h*ce + u;
dc  = n*ce*x - mu*c;
an = vpasolve(dx==0,dce==0,dc==0);
[an.x, an.c, an.ce]


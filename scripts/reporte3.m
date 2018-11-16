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
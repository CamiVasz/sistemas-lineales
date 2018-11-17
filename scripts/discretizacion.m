%Discretizacion
Ts = 0.5;

so = 0.3356; ro = 0.0683 ; uo= 0.12;
[A,B,C,D] = linmod('resistenciaBacteriasSal1_entrada',[so,ro],uo);
sys = ss(A,B,C,D);
sysD = c2d(sys,Ts);


%  step(sys(1))
%
%sys = sim('resistB_Compara','SaveState','on','StateSaveName','xout','SaveOutput','on','OutputSaveName','yout');

%--------------------------------------------------------
% Asignacion de polas

%[B,A] = numden(tf(sysD));
syms q0 q1 q2 r1 z 
Bp = (-0.003178*z - 0.002949);
Ap = (z^2 -1.794*z + 0.799);
Qc = q0*z^2 + q1*z + q2;
Pc = (z-r1)*(z-1);

PolosDes = [0.5 0.2 -.1 -.4];
Ad = (z-PolosDes(1))*(z-PolosDes(2))*(z-PolosDes(3))*(z-PolosDes(4));
vecAsig= coeffs(Ap*Pc + Bp*Qc - Ad,z);
solucionParam = solve(vecAsig==[0 0 0 0],q0,q1,q2,r1)



%-----------------------------------------------------
%Matriz de controlabilidad

Mc = [B A*B];
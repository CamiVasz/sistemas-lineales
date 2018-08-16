syms x y z
eq1 = y;
eq2 = x*z^2 - 4*x;
eq3 = x^2-1;
ans=vpasolve(eq1==0,eq2==0,eq3==0);
[ans.x,ans.y,ans.z]
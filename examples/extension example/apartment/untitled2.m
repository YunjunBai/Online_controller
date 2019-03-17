syms x y x1 x2  x3 x4 a b c d e f tau1 tau2 m1 m2 l1 l2 g
a=m2*l2^2;
b=m2*l1*l2*cos(x3);
c=(m1+m2)*l1^2;
d=m2*l1*l2*sin(x3);
e=m2*l2*g*sin(x1+x3);
f=(m1+m2)*l1*g*sin(x1);

eqn1 = a*(x+y)+b*(2*x+y)+c*x-d*x2^2-2*d*x2*x4+e+f==tau1;
eqn2 = b*x+d*x2^2+e+a*(x+y)==tau2;

sol = solve([eqn1, eqn2], [x, y]);
xSol = sol.x
ySol = sol.y


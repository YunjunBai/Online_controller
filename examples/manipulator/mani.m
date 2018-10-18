lb=[-pi/2,-pi/4,-pi/2,-pi/4];
ub=[pi/2,pi/4,pi/2,pi/4];
x0=[-pi/2,-pi/4,-pi/2,-pi/4];
[x1,mf]=fmincon(@func,x0,[],[],[],[],lb,ub);    

max=-mf
x1


% function f=func(x)
%  r_m=0.08;
% f=-((196*cos(x(1) + x(3)) - 98*sin(x(1) + x(3))*sin(x(3)) + 20*r_m*x(2)^2*cos(x(3)) + 196*sin(x(1))*sin(x(3)) + 15*r_m*x(2)^2*cos(x(3))^2 - 15*r_m*x(2)^2*sin(x(3))^2  + 98*cos(x(1) + x(3))*cos(x(3)) + 10*r_m*x(2)*r_m*x(4)*cos(x(3)) + 10*r_m*x(2)*r_m*x(4)*cos(x(3))^2 - 10*r_m*x(2)*r_m*x(4)*sin(x(3))^2)/(5*cos(x(3))^2 - 10) + (10*cos(x(3))*sin(x(3))*( 196*sin(x(1) + x(3)) - 196*sin(x(1)) + 20*r_m*x(2)^2*sin(x(3)) - 196*cos(x(3))*sin(x(1))  + 98*sin(x(1) + x(3))*cos(x(3)) + 10*r_m*x(2)*r_m*x(4)*sin(x(3)) + 15*r_m*x(2)^2*cos(x(3))*sin(x(3)) + 10*r_m*x(2)*r_m*x(4)*cos(x(3))*sin(x(3))))/(5*cos(x(3))^2 - 10)^2);
% 
% end

function f_1=func(x)
a1=0.25;
a2=100;
a3=10;
a4=0.1;
a5=100;
a6=0.01;
r_m=0.08;
m1=1;
m2=1;
l1=0.5;
l2=0.5;
tau_f1=a1*(tanh(a2*r_m*x(2))-tanh(a3*r_m*x(2)))+a4*tanh(a5*r_m*x(2))+a6*r_m*x(2);
tau_f2=a1*(tanh(a2*r_m*x(4))-tanh(a3*r_m*x(4)))+a4*tanh(a5*r_m*x(4))+a6*r_m*x(4);
f_1=-((l2*(tau_f2-tau_f1)+l1*tau_f2*cos(x(3)))/(l1^2*l2*m1+l1^2*l2*m2-l1^2*l2*m2*cos(x(3)^2)));
f_2=(-((l1^2*m1+l1^2*m2+l2^2*m2)*tau_f2 -l2^2*m2*tau_f1 + (2*l1*l2*m2*tau_f2-l1*l2*m2*tau_f1)*cos(x(3)))/(-0.0625*(cos(x(3))^2 )+ 0.125));
end
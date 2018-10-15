% syms x y x(1) x(2)  x(3) x(4)  tau1 tau2  tau_f1 tau_f2 
% m1=1;
% m2=1;
% l1=0.5;
% l2=0.5;
% g=9.8;
% a=m2*l2^2;
% b=m2*l1*l2*cos(x(3));
% c=(m1+m2)*l1^2;
% d=m2*l1*l2*sin(x(3));
% e=m2*l2*g*sin(x(1)+x(3));
% f=(m1+m2)*l1*g*sin(x(1));
% 
% 
% eqn1 = a*(x+y)+b*(2*x+y)+c*x-d*x(2)^2-2*d*x(2)*x(4)+e+f==tau1-tau_f1;
% eqn2 = b*x+d*x(2)^2+e+a*(x+y)==tau2-tau_f2;
% 
% sol = solve([eqn1, eqn2], [x, y]);
% xSol = sol.x
% ySol = sol.y

% % xSol =-(20*tau1 - 20*tau2 - 20*tau_f1 + 20*tau_f2 - 196*sin(x(1)) + 10*x(2)^2*sin(x(3)) - 20*tau2*cos(x(3)) + 20*tau_f2*cos(x(3)) + 98*sin(x(1) + x(3))*cos(x(3)) + 10*x(2)*x(4)*sin(x(3)) + 5*x(2)^2*cos(x(3))*sin(x(3)))/(5*(cos(x(3))^2 - 2))
% %  
% %  
% % ySol =(20*tau1 - 60*tau2 - 20*tau_f1 + 60*tau_f2 + 196*sin(x(1) + x(3)) - 196*sin(x(1)) + 20*x(2)^2*sin(x(3)) - 196*cos(x(3))*sin(x(1)) + 20*tau1*cos(x(3)) - 40*tau2*cos(x(3)) - 20*tau_f1*cos(x(3)) + 40*tau_f2*cos(x(3)) + 98*sin(x(1) + x(3))*cos(x(3)) + 10*x(2)*x(4)*sin(x(3)) + 15*x(2)^2*cos(x(3))*sin(x(3)) + 10*x(2)*x(4)*cos(x(3))*sin(x(3)))/(5*(cos(x(3))^2 - 2))
% %  
% % d1=diff(xSol, x(1))
% % d2=diff(xSol, x(2))
% % d3=diff(xSol, x(3))
% % d4=diff(xSol, x(4))
% % d5=diff(ySol, x(1))
% % d6=diff(ySol, x(2))
% % d7=diff(ySol, x(3))
% % d8=diff(ySol, x(4))
% % d1 =
% %  
% % (196*cos(x(1)) - 98*cos(x(1) + x(3))*cos(x(3)))/(5*cos(x(3))^2 - 10)
% %  
% %  
% % d2 =
% %  
% % -(20*x(2)*sin(x(3)) + 10*x(4)*sin(x(3)) + 10*x(2)*cos(x(3))*sin(x(3)))/(5*cos(x(3))^2 - 10)
% %  
% %  
% % d3 =
% %  
% % - (10*x(2)^2*cos(x(3)) - 98*sin(x(1) + x(3))*sin(x(3)) + 5*x(2)^2*cos(x(3))^2 - 5*x(2)^2*sin(x(3))^2 + 20*tau2*sin(x(3)) - 20*tau_f2*sin(x(3)) + 98*cos(x(1) + x(3))*cos(x(3)) + 10*x(2)*x(4)*cos(x(3)))/(5*cos(x(3))^2 - 10) - (10*cos(x(3))*sin(x(3))*(20*tau1 - 20*tau2 - 20*tau_f1 + 20*tau_f2 - 196*sin(x(1)) + 10*x(2)^2*sin(x(3)) - 20*tau2*cos(x(3)) + 20*tau_f2*cos(x(3)) + 98*sin(x(1) + x(3))*cos(x(3)) + 10*x(2)*x(4)*sin(x(3)) + 5*x(2)^2*cos(x(3))*sin(x(3))))/(5*cos(x(3))^2 - 10)^2
% %  
% %  
% % d4 =
% %  
% % -(10*x(2)*sin(x(3)))/(5*cos(x(3))^2 - 10)
% %  
% %  
% % d5 =
% %  
% % (196*cos(x(1) + x(3)) - 196*cos(x(1)) - 196*cos(x(1))*cos(x(3)) + 98*cos(x(1) + x(3))*cos(x(3)))/(5*cos(x(3))^2 - 10)
% %  
% %  
% % d6 =
% %  
% % (40*x(2)*sin(x(3)) + 10*x(4)*sin(x(3)) + 30*x(2)*cos(x(3))*sin(x(3)) + 10*x(4)*cos(x(3))*sin(x(3)))/(5*cos(x(3))^2 - 10)
% %  
% %  
% % d7 =
% %  
% % (196*cos(x(1) + x(3)) - 98*sin(x(1) + x(3))*sin(x(3)) + 20*x(2)^2*cos(x(3)) + 196*sin(x(1))*sin(x(3)) + 15*x(2)^2*cos(x(3))^2 - 15*x(2)^2*sin(x(3))^2 - 20*tau1*sin(x(3)) + 40*tau2*sin(x(3)) + 20*tau_f1*sin(x(3)) - 40*tau_f2*sin(x(3)) + 98*cos(x(1) + x(3))*cos(x(3)) + 10*x(2)*x(4)*cos(x(3)) + 10*x(2)*x(4)*cos(x(3))^2 - 10*x(2)*x(4)*sin(x(3))^2)/(5*cos(x(3))^2 - 10) + (10*cos(x(3))*sin(x(3))*(20*tau1 - 60*tau2 - 20*tau_f1 + 60*tau_f2 + 196*sin(x(1) + x(3)) - 196*sin(x(1)) + 20*x(2)^2*sin(x(3)) - 196*cos(x(3))*sin(x(1)) + 20*tau1*cos(x(3)) - 40*tau2*cos(x(3)) - 20*tau_f1*cos(x(3)) + 40*tau_f2*cos(x(3)) + 98*sin(x(1) + x(3))*cos(x(3)) + 10*x(2)*x(4)*sin(x(3)) + 15*x(2)^2*cos(x(3))*sin(x(3)) + 10*x(2)*x(4)*cos(x(3))*sin(x(3))))/(5*cos(x(3))^2 - 10)^2
% %  
% %  
% % d8 =
% %  
% % (10*x(2)*sin(x(3)) + 10*x(2)*cos(x(3))*sin(x(3)))/(5*cos(x(3))^2 - 10)
%  
% 
% 
%  
lb=[-pi/2,-pi/90,-pi/2,-pi/90];
ub=[pi/2,pi/90,pi/2,pi/90];
x0=[-pi/2,-1/pi,-pi/2,-1/pi];
[x1,mf]=fmincon(@func,x0,[],[],[],[],lb,ub);    

max=-mf
x1


% function f=func(x)
% 
% f=-((10*x(2)*sin(x(3)) + 10*x(2)*cos(x(3))*sin(x(3)))/(5*cos(x(3))^2 - 10));
% 
% end

function f_2=func(x)
a1=0.25;
a2=100;
a3=10;
a4=0.1;
a5=100;
a6=0.01;
m1=1;
m2=1;
l1=0.5;
l2=0.5;
tau_f1=a1*(tanh(a2*x(2))-tanh(a3*x(2)))+a4*tanh(a5*x(2))+a6*x(2);
tau_f2=a1*(tanh(a2*x(4))-tanh(a3*x(4)))+a4*tanh(a5*x(4))+a6*x(4);
f_1=-((l2*(tau_f2-tau_f1)+l1*tau_f2*cos(x(3)))/(l1^2*l2*m1+l1^2*l2*m2-l1^2*l2*m2*cos(x(3)^2)));
f_2=-(-((l1^2*m1+l1^2*m2+l2^2*m2)*tau_f2 -l2^2*m2*tau_f1 + (2*l1*l2*m2*tau_f2-l1*l2*m2*tau_f1)*cos(x(3)))/(-l1^2*l2^2*m2^2*cos(x(3)^2+l1^2*l2^2*m2^2+m1*l1^2*l2^2*m2)));
end
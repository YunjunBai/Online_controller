% % define an array of values for x y and z
% x2= -2:0.1:2; 
% x3 = -pi/2:pi/40:pi/2; 
% x4 = -2:0.1:2; 
% 
% % create a meshgrid
% [Xmesh,Ymesh,Zmesh] = meshgrid(x2,x3,x4);


a1=0.25;
a2=-100;
a3=10;
a4=0.1;
a5=100;
a6=0.01;
m1=1;
m2=1;
l1=0.5;
l2=0.5;
tau_f1=a1*(tanh(a2*x2)-tanh(a3*x2))+a4*tanh(a5*x2)+a6*x2;
tau_f2=a1*(tanh(a2*x4)-tanh(a3*x4))+a4*tanh(a5*x4)+a6*x4;
f_1=@(x2,x3,x4) (l2*(tau_f2-tau_f1)+l1*tau_f2*cos(x3))/(l1^2*l2*m1+l1^2*l2*m2-l1^2*l2*m2*cos(x3)^2);
f_2=@(x2, x3, x4) (-(l1^2*m1+l1^2*m2+l2^2*m2)*tau_f2 -l2^2*m2*tau_f1 + (2*l1*l2*m2*tau_f2-l1*l2*m2*tau_f1)*cos(x3))/(-l1^2*l2^2*m2^2*cos(x3^2+l1^2*l2^2*m2^2+m1*l1^2*l2^2*m2));

% launch sliceomatic
fsurf(f_1,[-2 2 -pi/2 pi/2 -2 2])



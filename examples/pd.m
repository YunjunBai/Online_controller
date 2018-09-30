% function pd
% syms u0 u1 x y z v;
% alphax=0.0204;
% alphay	=0.0242;
% betax=	0.0076;
% betay	=0.0168;
% k1=	0.0;
% k2	=2.0;
% k3	=0.8;
% k4	=0.5;
% m1	=0.00005;
% z0	=30.0;
% t	=12.5;
% r1=	15.0;
% r0=	10.0;
% c1=	0.0;
% c2=	0.0;
% c3=	0.0;
% d0=	1.0;
% 
% Gx	=((alphax * (k1 + ((1 - k1) * (z / (z + k2))))) - (betax * (k3 + (1 - k3) * (z / (z + k4)))));
% Gy=	((alphay * (1 - (d0 * (z / z0)))) - betay);
% Mxy=	(m1 * (1 - (z / z0)));
% scale	=50.0;
% 
% 
% 
% 
% f_x= scale * (((Gx - Mxy) * x) + c1 * x);
% f_y= scale * (((Mxy * x) + Gy * y ) + c2 * y);
% f_z = scale * (((z0 - z) / t) + c3 * z)*u0+scale * (((0 - z) / t) + c3 * z)*u1;
% f_v = scale * ((((Gx - Mxy) * x) + c1 * x) + (((Mxy * x) + Gy * y ) + c2 * y));
% 
% dxfx=diff(f_x,x)
% dyfx=diff(f_x,y)
% dzfx=diff(f_x,z)
% dvfx=diff(f_x,v)
% 
% dxfy=diff(f_y,x)
% dyfy=diff(f_y,y)
% dzfy=diff(f_y,z)
% dvfy=diff(f_y,v)
% 
% dxfz=diff(f_z,x)
% dyfz=diff(f_z,y)
% dzfz=diff(f_z,z)
% dvfz=diff(f_z,v)
% 
% dxfv=diff(f_v,x)
% dyfv=diff(f_v,y)
% dzfv=diff(f_v,z)
% dvfv=diff(f_v,v)
% 
% 
% end
% 
% dxfx =
%  
% z/12000 + (51*z)/(50*(z + 2)) - (19*z)/(250*(z + 1/2)) - 613/2000
%  
%  
% dzfx =
%  
% 50*x*((19*z)/(12500*(z + 1/2)^2) - (51*z)/(2500*(z + 2)^2) + 51/(2500*(z + 2)) - 19/(12500*(z + 1/2)) + 1/600000)
%  
%  
% dxfy =
%  
% 1/400 - z/12000
%  
%  
% dyfy =
%  
% 37/100 - (121*z)/3000
%  
%  
% dzfy =
%  
% - x/12000 - (121*y)/3000
%  
%  
%  
% dzfz =
%  
% - 4*u0 - 4*u1
%  
%  
%  
% dxfv =
%  
% (51*z)/(50*(z + 2)) - (19*z)/(250*(z + 1/2)) - 38/125
%  
%  
% dyfv =
%  
% 37/100 - (121*z)/3000
%  
%  
% dzfv =
%  
% 50*x*((19*z)/(12500*(z + 1/2)^2) - (51*z)/(2500*(z + 2)^2) + 51/(2500*(z + 2)) - 19/(12500*(z + 1/2)) + 1/600000) - (121*y)/3000 - x/12000
%  


lb=[0,0,0,0];
ub=[30,30,30,30];
x0=[0,0,0,0];
[x1,mf]=fmincon(@func,x0,[],[],[],[],lb,ub);    
max=-mf

x1


function f=func(x)
f=-(50*x(1)*((19*x(3))/(12500*(x(3) + 1/2)^2) - (51*x(3))/(2500*(x(3) + 2)^2) + 51/(2500*(x(3) + 2)) - 19/(12500*(x(3) + 1/2)) + 1/600000) - (121*x(2))/3000 - x(1)/12000);
end


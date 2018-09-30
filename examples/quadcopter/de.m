
% syms u w v p q r 
% 
% f_1=cos(b)*cos(c)*u+(-cos(a)*sin(c)+sin(a)*sin(b)*cos(c))*v+(sin(a)*sin(c)+cos(a)*sin(b)*cos(c))*w;
% f_2=cos(b)*sin(c)*u+(cos(a)*cos(c)+sin(a)*sin(b)*sin(c))*v+(-sin(a)*cos(c)+cos(a)*sin(b)*sin(c))*w;
% f_3=-sin(b)*u+sin(a)*cos(b)*v+cos(a)*cos(b)*w;
% f_4=p+(q*sin(a)+r*cos(a))*tan(b);
% f_5=q*cos(a)-r*sin(a);
% f_6=(q*sin(a)+r*cos(a))*sec(b);
% 
% p_14=diff(f_1,a);
% p_15=diff(f_1,b);
% p_16=diff(f_1,c);
% 
% p_24=diff(f_2,a);
% p_25=diff(f_2,b);
% p_26=diff(f_2,c);
% 
% p_34=diff(f_3,a);
% p_35=diff(f_3,b);
% p_36=diff(f_3,c);
% 
% p_44=diff(f_4,a);
% p_45=diff(f_4,b);
% p_46=diff(f_4,c);
% 
% p_55=diff(f_5,a);
% 
% p_64=diff(f_6,a);
% p_65=diff(f_6,b);
% 
% p_14,p_15,p_16,p_24,p_25,p_26,p_34,p_35,p_36,p_44,p_45,p_46,p_55,p_64,p_65
% 
% f=sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3))*sin(x(2));
% f=cos(x(1))*sin(x(3))-cos(x(3))*sin(x(1))*sin(x(2));
% f=cos(x(1))*cos(x(2))*cos(x(3));
% -f=cos(x(3))*sin(x(2));
% f=cos(x(2))*cos(x(2))*sin(x(1));
% f=cos(x(3))*sin(x(1))-cos(x(1))*sin(x(2))*sin(x(3));
% -f=cos(x(1))*cos(x(3))+sin(x(1))*sin(x(2))*sin(x(3));
% -f=cos(x(2))*cos(x(3));
% f=cos(x(3))*sin(x(1))-cos(x(1))*sin(x(2))*sin(x(3));
% f=cos(x(1))*cos(x(3))+sin(x(1))*sin(x(2))*sin(x(3));
% f=cos(x(1))*cos(x(2))*sin(x(3));
% -f=sin(x(x(2))*sin(x(3));
% f=cos(x(2))*sin(x(1))*sin(x(3));
% f=sin(x(1))*sin(x(3))+cos(x(1))*cos(x(2))*sin(x(2))
% -f=cos(x(1))*sin(x(3))-cos(x(3))*sin(x(1))*sin(x(2));
% f=cos(x(2))*cos(x(3));
% f=cos(x(1))*cos(x(2));
% -f=cos(x(2))*sin(x(1));
% f=cos(x(2))
% f=cos(x(1))*sin(x(2));
% f=sin(x(1))*sin(x(2));
% f=tan(x(2))*cos(x(1));
% -f=tan(x(2))*sin(x(1));
% f=(tan(x(2))^2 +1)*cos(x(1));
% f=(tan(x(2))^2 +1)*sin(x(1));
% f=cos(x(1));
% f=sin(x(1));
% f=cos(x(1))/cos(x(2))
%f=sin(x(1))/cos(x(2))
%f=sin(x(2))*cos(x(1))/(cos(x(2))*cos(x(2)))



ub=[60,60,60];
lb=[-20,-20,-20];
x0=[-20,-20,-20];
[x1,mf]=fmincon(@func,x0,[],[],[],[],lb,ub);
max=-mf

x1


function f=func(x)
f=-(sin(x(2))*sin(x(1))/(cos(x(2))*cos(x(2))));
end


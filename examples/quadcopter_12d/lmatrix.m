 syms u_0 u_1 u_2 u_3 y_3 y_4 y_5 y_6 y_7 y_8 y_9 y_10 y_11


     m=0.468;
     g=9.81;
     b=2.980*1e-6;
     l=0.225;
     d=1.140*1e-7;
     I_x=4.856*1e-3;
     I_y=4.856*1e-3;
     I_z=8.801*1e-3;
     f_t=b*(u_0*u_0+u_1*u_1+u_2*u_2+u_3*u_3);
     tau_x=b*l*(u_2*u_2-u_0*u_0);
     tau_y=b*l*(u_3*u_3-u_1*u_1);
     tau_z=d*(u_3*u_3+u_1*u_1-u_2*u_2-u_0*u_0);

    yy_0 = y_8*(sin(y_5)*sin(y_3)+cos(y_5)*cos(y_3)*sin(y_4)) - y_7*(cos(y_5)*sin(y_3)-cos(y_3)*sin(y_5)*sin(y_4)) +y_6*cos(y_3)*cos(y_4);
    yy_1 =-y_8*(cos(y_3)*sin(y_5)-cos(y_5)*sin(y_3)*sin(y_4)) + y_7*(cos(y_5)*cos(y_3)+sin(y_5)*sin(y_3)*sin(y_4)) +y_6*cos(y_4)*sin(y_3);
    yy_2 =y_8*cos(y_5)*cos(y_4)-y_6*sin(y_4)+y_7*cos(y_4)*sin(y_5);
    yy_3 =y_10*(sin(y_5)/cos(y_4))+y_11*(cos(y_5)/cos(y_4));
    yy_4 =y_10*cos(y_5)-y_11*sin(y_5);
    yy_5 =y_9+y_10*sin(y_5)*tan(y_4)+y_11*cos(y_5)*tan(y_4);
    yy_6 =-1/m *(sin(y_5)*sin(y_3)+cos(y_5)*cos(y_3)*sin(y_4)) ;
    yy_7 =-1/m *(cos(y_3)*sin(y_5)-cos(y_5)*sin(y_3)*sin(y_4)) ;
    yy_8 =-1/m *cos(y_5)*cos(y_4);
    yy_9 =(I_y-I_z)/I_x*y_10*y_11 + 1/I_x*tau_x;
    yy_10=(I_z-I_x)/I_y*y_9*y_11 + 1/I_y*tau_y;
    yy_11=(I_x-I_y)/I_z*y_9*y_10 +1/I_z *tau_z;


%  d03=diff(yy_0,y_3)
%  d04=diff(yy_0,y_4)
%  d05=diff(yy_0,y_5)
% d06=diff(yy_0,y_6)
% d07=diff(yy_0,y_7)
% d08=diff(yy_0,y_8)
   
% d24=diff(yy_2,y_4)
%  d25=diff(yy_2,y_5)
% d26=diff(yy_2,y_6)
% d27=diff(yy_2,y_7)
% d28=diff(yy_2,y_8)
% d34=diff(yy_3,y_4)
% d35=diff(yy_3,y_5)
% d310=diff(yy_3,y_10)
% d311=diff(yy_3,y_11)

% d45=diff(yy_4,y_5)
% d410=diff(yy_4,y_10)
% d411=diff(yy_4,y_11)
%  d54=diff(yy_5,y_4)
% d55=diff(yy_5,y_5)
% d59=diff(yy_5,y_9)
% d510=diff(yy_5,y_10)
% d511=diff(yy_5,y_11)
 d63=diff(yy_6,y_3)
 d64=diff(yy_6,y_4)
 d65=diff(yy_6,y_5)
 
  d73=diff(yy_7,y_3)
 d74=diff(yy_7,y_4)
 d75=diff(yy_7,y_5)

 d84=diff(yy_8,y_4)
 d85=diff(yy_8,y_5)
% d910=diff(yy_9,y_10)
% d911=diff(yy_9,y_11)
% d1010=diff(yy_10,y_9)
% d1011=diff(yy_10,y_11)
% 
% d1110=diff(yy_11,y_10)
% d1111=diff(yy_11,y_9)


function optDist
    % global bounds
%     fid = fopen('DistMatGlobal_1.txt','wt');
%     lb=[-pi/2,-pi/4,-pi/2,-pi/4];
%     ub=[pi/2,pi/4,pi/2,pi/4];
    % local changes
    lb=[-1.571	-0.471	0.257	-0.157	];
    ub=[1.885	0.6	1.712	0.571];
      fid = fopen('DistMatLocal_1.txt','wt');
    
    x0=[-pi/2,-pi/4,-pi/2,-pi/4];
    num = [10, 10, 10, 10];
    int = (ub-lb)./num;

    
    %[x1,mf]=fmincon(@func,x0,[],[],[],[],lb,ub);    

    % max=-mf
    % x1
%     a = zeros(1,4);
%     b = zeros(1,4);
    for i=1:num(1) % dim 1
        for j=1:num(2) % dim 2
            for k=1:num(3) % dim 3
                for l=1:num(4) % dim 4
                    a = [lb(1)+(i-1)*int(1),lb(2)+(j-1)*int(2),lb(3)+(k-1)*int(3),lb(4)+(l-1)*int(4)];
                    b = [a(1)+int(1),a(2)+int(2),a(3)+int(3),a(4)+int(4)];
                    if a(1)<lb(1) || a(2)<lb(2) || a(3)<lb(3) || a(4)<lb(4)
                        flag = 1;
                    end
                    if b(1)>ub(1) || b(2)>ub(2) || b(3)>ub(3) || b(4)>ub(4)
                        flag = 1;
                    end
                    [~,mf1]=fmincon(@func1,x0,[],[],[],[],a,b);
                    [~,mf2]=fmincon(@func2,x0,[],[],[],[],a,b);
                    max1 = -mf1;
                    max2 = -mf2;
                    write_to_file([0 abs(max1) 0 abs(max2) a b]);
                end
            end
        end
    end


    % function f=func(x)
    %  r_m=0.08;
    % f=-((196*cos(x(1) + x(3)) - 98*sin(x(1) + x(3))*sin(x(3)) + 20*r_m*x(2)^2*cos(x(3)) + 196*sin(x(1))*sin(x(3)) + 15*r_m*x(2)^2*cos(x(3))^2 - 15*r_m*x(2)^2*sin(x(3))^2  + 98*cos(x(1) + x(3))*cos(x(3)) + 10*r_m*x(2)*r_m*x(4)*cos(x(3)) + 10*r_m*x(2)*r_m*x(4)*cos(x(3))^2 - 10*r_m*x(2)*r_m*x(4)*sin(x(3))^2)/(5*cos(x(3))^2 - 10) + (10*cos(x(3))*sin(x(3))*( 196*sin(x(1) + x(3)) - 196*sin(x(1)) + 20*r_m*x(2)^2*sin(x(3)) - 196*cos(x(3))*sin(x(1))  + 98*sin(x(1) + x(3))*cos(x(3)) + 10*r_m*x(2)*r_m*x(4)*sin(x(3)) + 15*r_m*x(2)^2*cos(x(3))*sin(x(3)) + 10*r_m*x(2)*r_m*x(4)*cos(x(3))*sin(x(3))))/(5*cos(x(3))^2 - 10)^2);
    % 
    % end
    function write_to_file(x)
       fprintf(fid,'%0.3f\t',x);
       fprintf(fid,'\n');
    end
        

    function f_1=func1(x)
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

    function f_2=func2(x)
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

    fclose(fid);

end
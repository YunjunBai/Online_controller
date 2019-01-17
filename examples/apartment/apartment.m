% 
% function apartment
clear set
close all
figure
hold on

green = [0.7569    0.8667    0.7765];
		purple = [0.8196    0.6549    0.8471];
		orange = [0.9137    0.8275    0.3804];
		red = [0.9373    0.2980    0.2980];
		gray = [0.2471    0.2471    0.2471];

addpath(genpath('/var/tmp/online_algorithm/mfiles/'))
addpath(genpath('../../../mascot'))
load('data3');
%load('data_1');
% target set
lb=[9 0];
ub=lb+0.5;
% initial state
% x0=[0.6 0.6 0];
% 
% 
% % load controller from file
% %controller=StaticController('controller_global');
% 
% % simulate closed loop system
% y=x0;
v=[];
loop=3000;

% tau
tau=0.3;

% while(loop>0)
%     loop=loop-1;
%   
%   if (lb(1) <= y(end,1) & y(end,1) <= ub(1) &&...
%       lb(2) <= y(end,2) & y(end,2) <= ub(2))
%     break;
%   end 
% 
%   u=controller.control(y(end,:));
%   v=[v; u];
% 
%   [t x]=ode45(@unicycle_ode,[0 tau], y(end,:), odeset('abstol',1e-12,'reltol',1e-12),u);
%   
%   y=[y; x(end,:)];
% 
% end
% plot the vehicle domain
% colors
colors=get(groot,'DefaultAxesColorOrder');

% load the symbolic set containig obstacles
%obs=GridPoints('obstacles');
% obs=unique(obs(:,[1 2]),'rows');
% plot(obs(:,1),obs(:,2),'.');
% hold on

%target=GridPoints('target');
% target=unique(target(:,[1 2]),'rows');
% plot(target(:,1),target(:,2),'.','color',colors(2,:));
% hold on

% re_lazy=GridPoints('recomputation_lazy');
% re_lazy=unique(re_lazy(:,[1,2]),'rows');
% plot(re_lazy(:,1),re_lazy(:,2),'.','color',colors(2,:));




%winDomain_scots=StaticController('controller_global');
domain_scots=winDomain_scots.domain;
for k=1:size(domain_scots,1)
    x=domain_scots(k,1)-0.1;
    yy=domain_scots(k,2)-0.1;
    rectangle('Position',[x yy 0.2 0.2],'FaceColor',green,'EdgeColor',green');
%     rectangle('Position',[x y 0.2 0.2]);
end

%winDomain=StaticController('controller_scots2');
domain=winDomain.domain;
 for i=1:size(domain,1)
     x=domain(i,1)-0.1;
     yy=domain(i,2)-0.1;
     rectangle('Position',[x yy 0.2 0.2], 'FaceColor',purple,'EdgeColor',purple);
 end

% plot the real obstacles and target set
plot_domain

box on
axis([0 10 0 6])

lb=[9 0];
ub=lb+0.5;
% initial state
x0=[0.8 0.6 0];

% plot initial state  and trajectory

% for i=1:size(y,1)
    plot(y(:,1),y(:,2),'k.-')
    plot(y(1,1),y(1,2),'.','color',colors(5,:),'markersize',20)
% end


%set(gcf,'paperunits','centimeters','paperposition',[0 0 16 10],'papersize',[16 6])

% end

% **** Saving the figure in latex friendly way ******
set(gca,...
            'Units','inches',...
            'OuterPosition',[0 0 1.8 2],...
            'YTick',0:2:12,...
            'XTick',0:2:12,...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',9,...
            'FontName','Times')
        ylabel({'$x_2$'},...
                'FontUnits','points',...
                'interpreter','latex',...
                'FontSize',9,...
                'FontName','Times');%,...
%                 'Units','Normalized',...
%                 'Position',[-0.1,0.5,0])
        xlabel({'$x_1$'},...
                'FontUnits','points',...
                'interpreter','latex',...
                'FontSize',9,...
                'FontName','Times');%,...
%                 'Units','Normalized',...
%                 'Position',[0.5,-0.1,0])
    	
    	savefig('figure')
        
    print -depsc apartment.eps



function plot_domain
gray = [0.2471    0.2471    0.2471];
black = [0 0 0];
colors=get(groot,'DefaultAxesColorOrder');


alpha=0.2;

% v=[7.6 0;10 0;7.6 1.8; 10 1.8];
% patch('vertices',v,'faces',[1 2 4 3],'facea', alpha ,'edgec',colors(2,:));
% v=[4.2 0;7.4 0;4.2 1.8; 7.4 1.8]; % old disturbance region
v=[0.5 2;2.5 2;0.5 4; 2.5 4];
obj = patch('vertices',v,'faces',[1 2 4 3],'facea', alpha, 'facec',colors(3,:),'EdgeColor',black,'LineWidth',2);
hp = findobj(obj,'type','patch'); 
hatch(hp,45,black,'-',5,1)
hatch(hp,-45,black,'-',5,1)

v=[9 0; 9.5  0; 9 0.5; 9.5 .5];
obj = patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',colors(2,:),'EdgeColor',black,'LineWidth',2);
hp = findobj(obj,'type','patch'); 
hatch(hp,45,black,'-',5,1.5)

v=[4.0     0  ;4.2  0 ; 4.0  2 ; 4.2  2  ];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);
v=[5     1.8;7.6  1.8   ; 5   2.0    ; 7.6 2.0   ];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);
v=[4.0 3.4; 4.2 3.4; 4.0 6; 4.2 6];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);
v=[5 3.4; 7.6 3.4; 5 3.6; 7.6 3.6];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);
v=[7.4   0  ;7.6  0   ; 7.4   2.0   ; 7.6 2.0  ];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);
v=[7.4   3.4  ;7.6  3.4   ; 7.4   6.0    ; 7.6  6.0   ];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);
v=[8.4  1.8  ;10.0  1.8   ;8.4  2.0   ; 10.0  2.0  ];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);
v=[8.4   3.4  ;10.0    3.4   ; 8.4   3.6    ; 10.0   3.6   ];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',gray,'edgec',gray);

% Bounding box
rectangle('Position',[0 0 10 6],'EdgeColor',black,'LineWidth',2.5);

end

function dxdt = unicycle_ode(t,x,u)

  dxdt = zeros(3,1);
  c=atan(tan(u(2))/2);

  dxdt(1)=u(1)*cos(c+x(3))/cos(c);
  dxdt(2)=u(1)*sin(c+x(3))/cos(c);
  dxdt(3)=u(1)*tan(u(2));

end

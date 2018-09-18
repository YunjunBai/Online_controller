

function simple
clear set
close all


%% simulation

% target set
lb=[2.5 1.1];
ub=lb+0.5;
% initial state
x0=[0.6 0 -1.8];


% load controller from file
controller=StaticController('controller_1');

% simulate closed loop system
y=x0;
v=[];
loop=3000;

% tau
tau=0.3;

while(loop>0)
    loop=loop-1;
  
  if (lb(1) <= y(end,1) & y(end,1) <= ub(1) &&...
      lb(2) <= y(end,2) & y(end,2) <= ub(2))
    break;
  end 

  u=controller.control(y(end,:));
  v=[v; u];

  [t x]=ode45(@unicycle_ode,[0 tau], y(end,:), odeset('abstol',1e-12,'reltol',1e-12),u);
  
  y=[y; x(end,:)];

end

%% plot the vehicle domain
% colors
colors=get(groot,'DefaultAxesColorOrder');

% load the symbolic set containig obstacles
obs=GridPoints('obstacles');
obs=unique(obs(:,[1 2]),'rows');
plot(obs(:,1),obs(:,2),'.');
hold on

target=GridPoints('target');
target=unique(target(:,[1 2]),'rows');
plot(target(:,1),target(:,2),'.','color',colors(2,:));
hold on

% plot the real obstacles and target set
plot_domain

% plot initial state  and trajectory
 plot(y(:,1),y(:,2),'k.-')
 plot(y(1,1),y(1,2),'.','color',colors(5,:),'markersize',20)

box on
axis([0 3 0 3])

%set(gcf,'paperunits','centimeters','paperposition',[0 0 16 10],'papersize',[16 10])

end

function dxdt = unicycle_ode(t,x,u)

  dxdt = zeros(3,1);
  c=atan(tan(u(2))/2);

  dxdt(1)=u(1)*cos(c+x(3))/cos(c);
  dxdt(2)=u(1)*sin(c+x(3))/cos(c);
  dxdt(3)=u(1)*tan(u(2));

end

function plot_domain

colors=get(groot,'DefaultAxesColorOrder');

alpha=0.2;

v=[2.5 1.1; 3  1.1; 2.5 1.6; 3 1.6];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',colors(2,:),'edgec',colors(2,:));

v=[1.4    0 ;1.6  0   ; 1.4     2    ; 1.6  2   ];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'facec',colors(1,:),'edgec',colors(1,:));

v=[1.3 0.5; 2 0.5; 1.3 2.8; 2, 2.8];
patch('vertices',v,'faces',[1 2 4 3],'facea',alpha,'edgec',colors(3,:));
end

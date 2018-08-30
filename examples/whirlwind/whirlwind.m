%
% whirlwind.m
%
%
%
function whirlwind
clear set
close all


%% plot the  domain
% colors
colors=get(groot,'DefaultAxesColorOrder');

re_max=GridPoints('recomputation_max');
re_max=unique(re_max(:,[1 2]),'rows');
figure(1);
ploR2cartS(re_max(:,2),re_max(:,1));
%[re_max_x,re_max_y]=pol2cart(re_max(:,2),re_max(:,1));

plot(re_max_x,re_max_y,'.','color',colors(5,:));
%patch(re_max_x, re_max_y, 'white');
hold on
%plot_domain

box on
centers=[0,0];
B=0:0.5:10;
for i=1:length(B)
    viscircles(centers,B(i),'Color','black','LineStyle',':','LineWidth',0.01);
end  

% plot the "new disturbance region"

re_lazy=GridPoints('recomputation_lazy');
re_lazy=unique(re_lazy(:,[1 2]),'rows');
figure(2);
[re_lazy_x,re_lazy_y]=pol2cart(re_lazy(:,2),re_lazy(:,1));
plot(re_lazy_x,re_lazy_y,'.','color',colors(4,:));
%patch(re_lazy_x, re_lazy_y, 'green');
hold on
box on
centers=[0,0];
B=0:1:10;
for i=1:length(B)
    viscircles(centers,B(i),'Color','black','LineStyle',':','LineWidth',0.01);
end 

% plot the "new disturbance region"
hold on
end



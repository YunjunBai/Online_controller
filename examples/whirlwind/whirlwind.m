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
[re_max_x,re_max_y]=pol2cart(re_max(:,2),re_max(:,1));
plot(re_max_x,re_max_y,'.','color',colors(5,:));
hold on

re_lazy=GridPoints('recomputation_lazy');
re_lazy=unique(re_lazy(:,[1 2]),'rows');
[re_lazy_x,re_lazy_y]=pol2cart(re_lazy(:,2),re_lazy(:,1));
plot(re_lazy_x,re_lazy_y,'.','color',colors(4,:));
hold on

%plot_domain

box on
centers=[0,0];
radii=10.5;
viscircles(centers,radii);


% plot the "new disturbance region"



end



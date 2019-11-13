clc;clear;close all;tic;

s = 200; r = 0.5; t = 0.3;
al = 1:999;
al = al.*0.001;
al(1000)=0.9999;
dft_dal = zeros(1000,1);
dfs_dal = zeros(1000,1);
for i=1:1000
    dft_dal(i) = (s*t^2)/(al(i)*r - 1) - (al(i)*r*s*t^2)/(al(i)*r - 1)^2;
    dfs_dal(i) = (s*t)/(al(i)*r - 1) - (r*s*t*(al(i) - 1))/(al(i)*r - 1)^2;
end
al = al.';

figure
hold on
t = plot(al,dft_dal,'-');
s = plot(al,dfs_dal,'-');
hold off
ylabel('\partial F / \partial\alpha (W/m^2)','FontSize',24);xlabel('albedo (%)','FontSize',24);
title('Change in radiative flux w.r.t. albedo','FontSize',24);
legend([t,s],{'TOA','SFC'},'FontSize',20,'Location','northeast')
set(gca,'FontSize',20)

dal = zeros(size(al));
for i=1:size(al)
    dal(i)=al(i)-mean(al);
end

dft = dft_dal.*dal;
dfs = dfs_dal.*dal;

figure
hold on
t = plot(al,dft,'-');
s = plot(al,dfs,'-');
hold off
ylabel('\partial F (W/m^2)','FontSize',24);xlabel('albedo (%)','FontSize',24);
title('Change in radiative flux w.r.t. albedo','FontSize',24);
legend([t,s],{'TOA','SFC'},'FontSize',20,'Location','northeast')
set(gca,'FontSize',20)
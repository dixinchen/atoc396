clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-4/';

% lon = 205;
% lat = 70;
lon = 195;
lat = 70;
W = cos(lat'*pi/180);

% coord index for 2.5*2.5 resolution
lon_k = ncread([result_path,'Feedback_kernel_cld.nc'],'longitude');
lat_k = ncread([result_path,'Feedback_kernel_cld.nc'],'latitude');
idx_lon_k = find(lon_k==lon);
idx_lat_k = find(lat_k==lat);

% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'latitude');
idx_lon_nn = find(lon_nn==lon);
idx_lat_nn = find(lat_nn==lat);

trend = "deseasonal";
% trend = "interannual";

% var = "wv";
var = "fal";

month = 9;
num = 10;
%%
clim_mean = ncread([gen_path, 'data_1_m/Test2007_2016/fal_mean_10.nc'],'fal_mean');
% clim_mean = ncread([gen_path, 'data_1_m/Test2007_2016/tcwv_mean_10.nc'],'tcwv_mean');
clim_mean = squeeze(clim_mean(idx_lon_nn, idx_lat_nn, month));

x = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
% x = ncread([gen_path,'data_1_m/GTCWV_1_m38.nc'],'tcwv');
if trend == "deseasonal"
%     x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
%     x = squeeze(x(idx_lon_nn, idx_lat_nn, month:12:12*37+month));
    x = squeeze(x(idx_lon_nn, idx_lat_nn, :));
end
if trend == "interannual"
%     x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
%     x = squeeze(x(idx_lon_nn, idx_lat_nn, month:12:12*37+month));
    x = squeeze(x(idx_lon_nn, idx_lat_nn, :));
    x = x.';
    X = 1:456;
    p = polyfit(X',x',1);
    x_dt = x'-p(1)*X'-p(2);
    x = x_dt;
end

dr = ncread([gen_path,'data_1_m/GSSRC_1_m38.nc'],'ssrc');
if trend == "deseasonal"
%     ssr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
%     ssr = squeeze(dr(idx_lon_nn, idx_lat_nn, month:12:12*37+month));
    ssr = squeeze(dr(idx_lon_nn, idx_lat_nn, :))./2.628e+6;
end
if trend == "interannual"
%     dr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
%     dr = squeeze(dr(idx_lon_nn, idx_lat_nn, month:12:12*37+month));
    dr = squeeze(dr(idx_lon_nn, idx_lat_nn, :))./2.628e+6;
    dr = dr.';
    X = 1:456;
    p = polyfit(X',dr',1);
    ssr = dr'-p(1)*X'-p(2);
end



figure;
plot(x,ssr,'.');
ylabel('SSRC (W/m^2)','FontSize',24);xlabel('fal','FontSize',24);
title('38*12 monthly-mean at (70, 195), '+trend,'FontSize',24);
% ylabel('dRw (W/m^2)','FontSize',24);xlabel('tcwv','FontSize',24);
% title('2007-2016, Sept, SW F^-_{sfc,w}/tcwv, '+trend,'FontSize',24);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',20)

% a = x(1);
% b = ssr(1);
% 
% for i=2:size(dr)
%     if (ssr(i) == min(ssr)) || (x(i) == min(x)) || (ssr(i)<(3*10^6)) || (x(i)<0.1)
%         continue
%     end
%     a(end+1) = x(i);
%     b(end+1) = ssr(i);
% end
% 
% plot(a,b,'+');
% hold off

%%
num = 20;
dr = ncread('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test1997_2016/dfal_1_m20.nc','d_fal');
if trend == "deseasonal"
    dfal = squeeze(dr(idx_lon_k, idx_lat_k,month,:));
end
if trend == "interannual"
    dr = squeeze(dr(idx_lon_k, idx_lat_k,month,:));
    dr = dr.';
    X = 1:num;
    p = polyfit(X',dr',1);
    dr_dt = dr'-p(1)*X'-p(2);
    dfal = dr_dt;
end



dr = ncread([result_path20,'Feedback_NN_cld_20yrs.nc'],'dRa');
if trend == "deseasonal"
    dr_nn = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
end
if trend == "interannual"
    dr = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
    dr = dr.';
    X = 1:num;
    p = polyfit(X',dr',1);
    dr_dt = dr'-p(1)*X'-p(2);
    dr_nn = dr_dt;
end

fit_x = linspace(max(dfal), min(dfal), 1000);
fit_dr_nn = polyval(polyfit(dfal, dr_nn, 2),fit_x);

figure;
hold on
plot(dfal,dr_nn,'.','MarkerSize',15);
plot(fit_x,fit_dr_nn, 'r--');
% if trend == "deseasonal"
%     mean = xline(clim_mean);
% end
hold off
ylabel('dRa (W/m^2)','FontSize',24);xlabel('dfal','FontSize',24);
title('2007-2016, Jun, SW F^-_{sfc,a}/fal, '+trend,'FontSize',24);
% ylabel('dRw (W/m^2)','FontSize',24);xlabel('dtcfal (kgm^{-2})','FontSize',24);
% title('2007-2016, Sept, SW F^-_{sfc,w}/tcwv, '+trend,'FontSize',24);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',20)
% if trend == "deseasonal"
%     legend([nn,mean],{'Kernel','NN','Clim. mean'},'FontSize',20);
% end
% if trend == "interannual"
%     legend([k,nn],{'Kernel','NN'},'FontSize',20);
% end
% saveas(gcf,figure_path+trend+'-drx-dx-'+var+'.png');

%%

clim_mean = ncread([gen_path, 'data_1_m/Test2007_2016/fal_mean_10.nc'],'fal_mean');
% clim_mean = ncread([gen_path, 'data_1_m/Test2007_2016/tcwv_mean_10.nc'],'tcwv_mean');
clim_mean = squeeze(clim_mean(idx_lon_nn, idx_lat_nn, month));

x = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
% x = ncread([gen_path,'data_1_m/GTCWV_1_m38.nc'],'tcwv');
if trend == "deseasonal"
    x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
end
if trend == "interannual"
    x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
    x = x.';
    X = 1:num;
    p = polyfit(X',x',1);
    x_dt = x'-p(1)*X'-p(2);
    x = x_dt;
end

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
if trend == "deseasonal"
    dr_k = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
end
if trend == "interannual"
    dr = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
    dr = dr.';
    X = 1:num;
    p = polyfit(X',dr',1);
    dr_dt = dr'-p(1)*X'-p(2);
    dr_k = dr_dt;
end

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
if trend == "deseasonal"
    dr_nn = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,month,:));
end
if trend == "interannual"
    dr = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,month,:));
    dr = dr.';
    X = 1:num;
    p = polyfit(X',dr',1);
    dr_dt = dr'-p(1)*X'-p(2);
    dr_nn = dr_dt;
end

fit_x = linspace(max(x), min(x), 1000);
fit_dr_k = polyval(polyfit(x, dr_k, 1),fit_x);
fit_dr_nn = polyval(polyfit(x, dr_nn, 2),fit_x);

figure;
hold on
k = plot(x,dr_k, '.','MarkerSize',15);
nn = plot(x,dr_nn,'.','MarkerSize',15);
plot(fit_x,fit_dr_k, 'b--');
plot(fit_x,fit_dr_nn, 'r--');
if trend == "deseasonal"
    mean = xline(clim_mean);
end
hold off
ylabel('dRa (W/m^2)','FontSize',24);xlabel('fal','FontSize',24);
title('2007-2016, Jun, SW F^-_{sfc,a}/fal, '+trend,'FontSize',24);
% ylabel('dRw (W/m^2)','FontSize',24);xlabel('tcwv (kgm^{-2})','FontSize',24);
% title('2007-2016, Sept, SW F^-_{sfc,w}/tcwv, '+trend,'FontSize',24);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',20)
if trend == "deseasonal"
    legend([k,nn,mean],{'Kernel','NN','Clim. mean'},'FontSize',20);
end
if trend == "interannual"
    legend([k,nn],{'Kernel','NN'},'FontSize',20);
end
saveas(gcf,figure_path+trend+'-drx-dx-'+var+'.png');
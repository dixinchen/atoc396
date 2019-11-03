clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';

lon = 205;
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

month = 9;
num = 10;

clim_mean = ncread([gen_path, 'data_1_m/Test2007_2016/fal_mean_10.nc'],'fal_mean');
clim_mean = squeeze(clim_mean(idx_lon_nn, idx_lat_nn, month));

x = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
dr = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
dra_k = dr;
% X = 1:num;
% p = polyfit(X',dr',1);
% dr_dt = dr'-p(1)*X'-p(2);
% dr = dr_dt.';

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
dr = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,month,:));
% X = 1:num;
% p = polyfit(X',dr',1);
% dr_dt = dr'-p(1)*X'-p(2);
% dr = dr_dt.';

fit_x = linspace(max(x), min(x), 1000);
fit_dr_k = polyval(polyfit(x, dr_k, 1),fit_x);
fit_dr_nn = polyval(polyfit(x, dr_nn, 2),fit_x);

figure;
hold on
k = plot(x,dr_k, '.','MarkerSize',15);
nn = plot(x,dr_nn,'.','MarkerSize',15);
truth = plot(x,dr_true,'k.','MarkerSize',15);
plot(fit_x,fit_dr_k, 'b--');
plot(fit_x,fit_dr_nn, 'r--');
mean = xline(clim_mean);
hold off
ylabel('dR','FontSize',16);xlabel('fal (%)','FontSize',16);title('2007-2016, Sept, SFC, SW, deseasonaled','FontSize',16)
legend([k,nn,truth,mean],{'Kernel dRa','NN dRa','EARi dSSR','Clim. mean'},'FontSize',16);


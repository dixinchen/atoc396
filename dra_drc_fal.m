clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';

m = ["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"];

trend = "deseasonal";
% trend = "interannual";
var = "fal";
month = 4;
num = 10;
lon = 75; lat = 75;
% lon = 75; lat = 70;

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

% get albedo data
data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
x = squeeze(data(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));

data = ncread([gen_path,'data_1_m/anomaly/dfal_1_m38.nc'],'d_fal');
dfal = squeeze(data(idx_lon_nn, idx_lat_nn, month, 29:end));

% get NN dra data (TOA SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRa');
dra_nn = squeeze(data(1,1,idx_lon_nn, idx_lat_nn, month, :));

% get NN drc data (TOA SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRc');
drc_nn = squeeze(data(1,1,idx_lon_nn, idx_lat_nn, month, :));


% get kernel dra data (TOA SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRa');
dra_k = squeeze(data(1,1,idx_lon_k, idx_lat_k, month, :));

% get kernel drc data (TOA SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRc');
drc_k = squeeze(data(1,1,idx_lon_k, idx_lat_k, month, :));

% get TSR data
data = ncread([gen_path,'data_1_m/anomaly/dtsr_1_m38.nc'],'d_tsr');
tsr = squeeze(data(idx_lon_nn, idx_lat_nn, month, 29:end));

% regression
fit_x = linspace(max(x), min(x), 1000);
fit_dra_nn = polyval(polyfit(x, dra_nn, 2),fit_x);
fit_drc_nn = polyval(polyfit(x, drc_nn, 2),fit_x);
fit_dra_k = polyval(polyfit(x, dra_k, 1),fit_x);
fit_drc_k = polyval(polyfit(x, drc_k, 1),fit_x);
fit_tsr = polyval(polyfit(x, tsr, 2),fit_x);

figure
hold on
a_k = plot(x,dra_k,'ro','MarkerSize',6);
c_k = plot(x,drc_k,'r^','MarkerSize',6);
a_nn = plot(x,dra_nn,'bo','MarkerSize',6);
c_nn = plot(x,drc_nn,'b^','MarkerSize',6);
plot(fit_x,fit_dra_nn,'r--','LineWidth',1);
plot(fit_x,fit_drc_nn,'r--','LineWidth',1);
plot(fit_x,fit_dra_k,'b--','LineWidth',1);
plot(fit_x,fit_drc_k,'b--','LineWidth',1);
true = plot(x,tsr,'k.','MarkerSize',14);
plot(fit_x,fit_tsr,'k-','LineWidth',1);
hold off
legend([a_k,c_k,a_nn,c_nn,true],{'dRa kernel','dRc kernel','dRa NN','dRc NN','dTSR (EARi)'},'FontSize',20)
ylabel("dR [W/m^2]",'FontSize',20);xlabel("fal [%]",'FontSize',20);
set(gca,'FontSize',20)


% regression
fit_x = linspace(max(dfal), min(dfal), 1000);
fit_dra_nn = polyval(polyfit(dfal, dra_nn, 2),fit_x);
fit_drc_nn = polyval(polyfit(dfal, drc_nn, 2),fit_x);
fit_dra_k = polyval(polyfit(dfal, dra_k, 1),fit_x);
fit_drc_k = polyval(polyfit(dfal, drc_k, 1),fit_x);
fit_tsr = polyval(polyfit(dfal, tsr, 2),fit_x);

figure
hold on
a_k = plot(dfal,dra_k,'ro','MarkerSize',6);
c_k = plot(dfal,drc_k,'r^','MarkerSize',6);
a_nn = plot(dfal,dra_nn,'bo','MarkerSize',6);
c_nn = plot(dfal,drc_nn,'b^','MarkerSize',6);
plot(fit_x,fit_dra_nn,'r--','LineWidth',1);
plot(fit_x,fit_drc_nn,'r--','LineWidth',1);
plot(fit_x,fit_dra_k,'b--','LineWidth',1);
plot(fit_x,fit_drc_k,'b--','LineWidth',1);
true = plot(dfal,tsr,'k.','MarkerSize',14);
plot(fit_x,fit_tsr,'k-','LineWidth',1);
hold off
grid on
legend([a_k,c_k,a_nn,c_nn,true],{'dRa kernel','dRc kernel','dRa NN','dRc NN','dTSR (EARi)'},'FontSize',20)
ylabel("dR [W/m^2]",'FontSize',20);xlabel("\Delta fal [%]",'FontSize',20);
set(gca,'FontSize',20)
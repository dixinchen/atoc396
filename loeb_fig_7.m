clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';

m = ["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"];

trend = "deseasonal";
% trend = "interannual";
var = "fal";
month = 6:8;
num = 10;
lon = 75; lat = 70;
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

% get NN dra data (TOA SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRa');
dra_nn = squeeze(data(2,1,idx_lon_nn, idx_lat_nn, month, :));
dra_nn = mean(dra_nn,1);

% get NN drc data (TOA SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRc');
drc_nn = squeeze(data(2,1,idx_lon_nn, idx_lat_nn, month, :));
drc_nn = mean(drc_nn,1);

% get kernel dra data (TOA SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRa');
dra_k = squeeze(data(1,1,idx_lon_k, idx_lat_k, month, :));
dra_k = mean(dra_k,1);

% get kernel drc data (TOA SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRc');
drc_k = squeeze(data(1,1,idx_lon_k, idx_lat_k, month, :));
drc_k = mean(drc_k,1);

total = drc_nn+dra_nn;

data = ncread([gen_path,'data_1_m/anomaly/dfal_1_m38.nc'],'d_fal');
dfal = squeeze(data(idx_lon_nn, idx_lat_nn, month, 29:end));
dfal = mean(dfal,1);

% get TSR data
data = ncread([gen_path,'data_1_m/anomaly/dtsr_1_m38.nc'],'d_tsr');
tsr = squeeze(data(idx_lon_nn, idx_lat_nn, month, 29:end));
tsr = mean(tsr,1);

x = linspace(2007, 2016, 10);

figure
hold on
a_nn = plot(x,dra_nn,'-^','MarkerSize',7,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','blue');
c_nn = plot(x,drc_nn,'-o','MarkerSize',7,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor','red');
% true = plot(x,dfal,'k-s','MarkerSize',7,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','black');
yline(0);
hold off
% legend([a_k,c_k,a_nn,c_nn,true],{'dRa kernel','dRc kernel','dRa NN','dRc NN','dTSR (EARi)'},'FontSize',20)
% ylabel("dR [W/m^2]",'FontSize',20);xlabel("fal [%]",'FontSize',20);
set(gca,'FontSize',20)
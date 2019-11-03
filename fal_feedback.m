clc;clear;close all;tic;

num = 12*10;
gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';

%% set month and coordinate of interest
month = 6;
lon = 205;
lat = 80;
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

%% load climatological mean
cli_mean = ncread([gen_path, 'data_1_m/Test2007_2016/fal_mean_10.nc'],'fal_mean');
cli_mean = cli_mean(idx_lon_nn, idx_lat_nn, :).*W;

%% load EARi data
fal = ncread([gen_path, 'data_1_m/GFAL_1_m38.nc'], 'fal');
fal = fal(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month);
fal = reshape(fal,[1 1 10]).*repmat(W,[1 1 10]);
fal = squeeze(nansum(fal,1));
fal = fal.';
% X = 1:10;
% p = polyfit(X',fal',1);
% fal_dt = fal'-p(1)*X'-p(2);

%% load true data

% TOA SW
dR = ncread([gen_path,'data_1_m/Test2007_2016/dtsr_1_m10.nc'],'d_tsr');
dR = squeeze(dR(idx_lon_nn, idx_lat_nn, month, :));
dtsr = reshape(dR,[1 1 10]).*repmat(W,[1 1 10]);
dtsr = squeeze(nansum(dtsr,1));
dtsr = dtsr.';

% SFC SW
dR = ncread([gen_path,'data_1_m/Test2007_2016/dssr_1_m10.nc'],'d_ssr');
dR = squeeze(dR(idx_lon_nn, idx_lat_nn, month, :));
dssr = reshape(dR,[1 1 10]).*repmat(W,[1 1 10]);
dssr = squeeze(nansum(dssr,1));
dssr = dssr.';

%% load kernel data
dR = ncread([result_path,'Feedback_kernel_cld.nc'],'dR');
dR = squeeze(dR(2, 1, idx_lon_k, idx_lat_k, month, :));
dR_k = reshape(dR,[1 1 10]).*repmat(W,[1 1 10]);
dR_k = squeeze(nansum(dR_k,1));
dR_k = dR_k.';

%% load NN data
dR = ncread([result_path,'Feedback_NN_cld.nc'],'dR');
dR = squeeze(dR(2, 1, idx_lon_nn, idx_lat_nn, month, :));
dR_nn = reshape(dR,[1 1 10]).*repmat(W,[1 1 10]);
dR_nn = squeeze(nansum(dR_nn,1));
dR_nn = dR_nn.';

%%
clearvars -except fal dtsr dR_k dR_nn cli_mean month

%% plot
figure;
plot(fal, dtsr, '.');
hold on
plot(fal, dR_k, '.');
hold on
plot(fal, dR_nn, '.');
hold on
xline(cli_mean(month));
hold on
legend('EARi', 'Kernel', 'NN', 'Clim. mean');
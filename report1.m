clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';

m = ["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"];

trend = "deseasonal";

var = "fal";

month = 6;
num = 10;

lon = 75; lat = 75;

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

%% 
% get SSR data
data = ncread([gen_path,'data_1_m/Test2007_2016/dssr_1_m10.nc'],'d_ssr');
ssr = squeeze(data(idx_lon_nn, idx_lat_nn, month,:));
ssr = ssr(:);

% get TSR data
data = ncread([gen_path,'data_1_m/Test2007_2016/dtsr_1_m10.nc'],'d_tsr');
tsr = squeeze(data(idx_lon_nn, idx_lat_nn, month,:));
tsr = tsr(:);

%%
% get NN method dra data

% % get dfal data
data = ncread([gen_path,'data_1_m/Test2007_2016/dfal_1_m10.nc'],'d_fal');
x = squeeze(data(idx_lon_nn, idx_lat_nn, month,:));
x = x(:);

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
dr_k = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
dr_k = dr_k(:);

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
dr_nn = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,month,:));
dr_nn = dr_nn(:);

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRc');
drc_k = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
drc_k = drc_k(:);

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRc');
drc_nn = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,month,:));
drc_nn = drc_nn(:);

fit_x = linspace(max(x), min(x), 1000);
fit_drk = polyval(polyfit(x, dr_k, 1),fit_x);
fit_drnn = polyval(polyfit(x, dr_nn, 2),fit_x);
fit_drck = polyval(polyfit(x, drc_k, 1),fit_x);
fit_drcnn = polyval(polyfit(x, drc_nn, 2),fit_x);

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
dr = squeeze(dr(1,1,idx_lon_k, idx_lat_k,month,:));
dr_kT = dr(:);

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
dr = squeeze(dr(1,1,idx_lon_nn, idx_lat_nn,month,:));
dr_nnT = dr(:);

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRc');
drc_kT = squeeze(dr(1,1,idx_lon_k, idx_lat_k,month,:));
drc_kT = drc_kT(:);

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRc');
drc_nnT = squeeze(dr(1,1,idx_lon_nn, idx_lat_nn,month,:));
drc_nnT = drc_nnT(:);

fit_x = linspace(max(x), min(x), 1000);
fit_drkT = polyval(polyfit(x, dr_kT, 1),fit_x);
fit_drnnT = polyval(polyfit(x, dr_nnT, 2),fit_x);
fit_drckT = polyval(polyfit(x, drc_kT, 1),fit_x);
fit_drcnnT = polyval(polyfit(x, drc_nnT, 2),fit_x);

figure
set(gca,'FontSize',20)
grid on

subplot(2,1,1)
hold on
plot(fit_x,fit_drnn,'--');
plot(fit_x,fit_drcnn,'--');
plot(x, dr_nn, 'o');
plot(x, drc_nn, '^');
hold off
% title("SFC",'FontSize',24);
ylabel('dR_x (W/m^2)','FontSize',24);xlabel('\Delta fal','FontSize',24);
% legend([ann,cnn],["Albedo","Cloud"],'FontSize',20)

% subplot(2,1,2)
% hold on
% plot(fit_x,fit_drnnT,'--');
% plot(fit_x,fit_drcnnT,'--');
% ann = plot(x, dr_nnT, 'o');
% cnn = plot(x, drc_nnT, '^');
% hold off
% ylabel('dR_x (W/m^2)','FontSize',24);xlabel('\Delta fal','FontSize',24);
% % title(strcat("Albedo feedback anomaly in ",num2str(m(month))),'FontSize',24);
% title("TOA",'FontSize',24);
% % legend([ann,cnn],{'Albedo','Cloud'},'FontSize',20)

% saveas(gcf,strcat(figure_path, "dRa_fal_",num2str(m(month)),".png"));
% saveas(gcf,strcat(figure_path, "dRa_dfal_AMJJ.png"));
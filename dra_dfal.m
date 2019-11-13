clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';

m = ["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"];

trend = "deseasonal";
% trend = "interannual";

var = "fal";

month = 3;
num = 10;

lon = 75; lat = 75;
% lon = 195; lat = 70;

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
% % get albedo data
data = ncread([gen_path,'data_1_m/Test2007_2016/dfal_1_m10.nc'],'d_fal');
% x = zeros(4,38);
% for i=1:38
%     for j = 4:7
%         x(j-3,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j));
%     end
% end 
x = squeeze(data(idx_lon_nn, idx_lat_nn, 4:7,:));
x = x(:);
% data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
% x = zeros(4,10);
% for i=1:10
%     for j = 4:7
%         x(j-3,i)=squeeze(data(idx_lon_nn, idx_lat_nn, 12*28+(i-1)*12+j));
%     end
% end 
% x = x(:);


% get SSR data
data = ncread([gen_path,'data_1_m/Test2007_2016/dssr_1_m10.nc'],'d_ssr');
% ssr = zeros(4,38);
% for i=1:38
%     for j = 4:7
%         ssr(j-3,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
%     end
% end
ssr = squeeze(data(idx_lon_nn, idx_lat_nn, 4:7,:));
ssr = ssr(:);

% get TSR data
data = ncread([gen_path,'data_1_m/Test2007_2016/dtsr_1_m10.nc'],'d_tsr');
% tsr = zeros(4,38);
% for i=1:38
%     for j = 4:7
%         tsr(j-3,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
%     end
% end
tsr = squeeze(data(idx_lon_nn, idx_lat_nn, 4:7,:));
tsr = tsr(:);

% % plot dTSR and dSSR against dFAL
% figure;
% hold on
% t = plot(x, tsr, '.');
% s = plot(x, ssr, '.');
% hold off
% ylabel('net solar radiation (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
% title('net SW flux anomaly over fal anomaly at (75, 75), Apr-Jul','FontSize',24);
% legend([t,s],{'TSR','SSR'},'FontSize',20)
% set(gca,'FontSize',20)
% set(gcf, 'Position', [350, 350, 650, 400]);
% % saveas(gcf,strcat(figure_path, "amjj.png"));

%%
% get NN method dra data

% % get GFAL data AMJJ for 10 yrs
% data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
% % x=squeeze(data(idx_lon_nn, idx_lat_nn, 12*28+month:12:end));
% % x = x(:);
% x = zeros(10,4);
% for j = 4:7
%     x(:,j-3) = squeeze(data(idx_lon_nn, idx_lat_nn, 12*28+j:12:12*37+j));
% end
% x = x(:);


month = 4:7;
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

fit_x = linspace(max(x), min(x), 1000);
fit_drk = polyval(polyfit(x, dr_k, 1),fit_x);
fit_drnn = polyval(polyfit(x, dr_nn, 2),fit_x);

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
dr = squeeze(dr(1,1,idx_lon_k, idx_lat_k,month,:));
dr_kT = dr(:);

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
dr = squeeze(dr(1,1,idx_lon_nn, idx_lat_nn,month,:));
dr_nnT = dr(:);

fit_x = linspace(max(x), min(x), 1000);
fit_drkT = polyval(polyfit(x, dr_kT, 1),fit_x);
fit_drnnT = polyval(polyfit(x, dr_nnT, 2),fit_x);

figure
hold on
plot(fit_x,fit_drk,'r--');
plot(fit_x,fit_drnn,'r--');
k = plot(x, dr_k, 'r^');
nn = plot(x, dr_nn, 'ro');
plot(fit_x,fit_drkT,'b--');
plot(fit_x,fit_drnnT,'b--');
kT = plot(x, dr_kT, 'b^');
nnT = plot(x, dr_nnT, 'bo');
hold off
ylabel('dRa (W/m^2)','FontSize',24);xlabel('dFAL','FontSize',24);
% title(strcat("Albedo feedback anomaly in ",num2str(m(month))),'FontSize',24);
title(strcat("Albedo feedback anomaly in AMJJ"),'FontSize',24);
legend([k,nn,kT,nnT],{'SFC Kernel','SFC NN','TOA Kernel','TOA NN'},'FontSize',20)
set(gca,'FontSize',20)
% saveas(gcf,strcat(figure_path, "dRa_fal_",num2str(m(month)),".png"));
% saveas(gcf,strcat(figure_path, "dRa_dfal_AMJJ.png"));
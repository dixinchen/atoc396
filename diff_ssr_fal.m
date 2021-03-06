clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';

m = ["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"];

trend = "deseasonal";
% trend = "interannual";

var = "fal";

month = 5:7;
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

% % get albedo data
% data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
% x = zeros(3,38);
% for i=1:38
%     for j = 5:7
%         x(j-4,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j));
%     end
% end 
% x = x(:);
% 
% 
% % get SSR data
% data = ncread([gen_path,'data_1_m/GSSR_1_m38.nc'],'ssr');
% ssr = zeros(3,38);
% for i=1:38
%     for j = 5:7
%         ssr(j-4,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
%     end
% end
% ssr = ssr(:);
% 
% % get TSR data
% data = ncread([gen_path,'data_1_m/GTSR_1_m38.nc'],'tsr');
% tsr = zeros(3,38);
% for i=1:38
%     for j = 5:7
%         tsr(j-4,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
%     end
% end
% tsr = tsr(:);

% get albedo data
data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
x = zeros(month(2)-month(1)+1,38);
for i=1:38
    for j = month
        x(j-(month(1)-1),i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j));
    end
end 
x = x(:);


% get SSR data
data = ncread([gen_path,'data_1_m/GSSR_1_m38.nc'],'ssr');
ssr = zeros(month(2)-month(1)+1,38);
for i=1:38
    for j = month
        ssr(j-(month(1)-1),i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
    end
end
ssr = ssr(:);

% get TSR data
data = ncread([gen_path,'data_1_m/GTSR_1_m38.nc'],'tsr');
tsr = zeros(month(2)-month(1)+1,38);
for i=1:38
    for j = month
        tsr(j-(month(1)-1),i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
    end
end
tsr = tsr(:);


%%
% fit TSR and SSR and get dTSR/dFAL and dSSR/dFAL
fit_x = linspace(max(x), min(x), 1000);
fit_ft = polyval(polyfit(x, tsr, 2),fit_x);
fit_fs = polyval(polyfit(x, ssr, 2),fit_x);

pt = polyfit(x, tsr, 2);
ps = polyfit(x, ssr, 2);

%%
syms ft fs al
ft = pt(1)*al^2+pt(2)*al+pt(3);
dft_eqn = diff(ft,al);
fs = ps(1)*al^2+ps(2)*al+ps(3);
dfs_eqn = diff(fs,al);

dft = zeros(size(tsr));
dfs = zeros(size(tsr));
for i=1:size(tsr)
%     dft(i) = 7289531910399253/35184372088832 - (4703651656153431*x(i))/4398046511104;
%     dfs(i) = - (7755799267377777*x(i))/17592186044416 - 6196388452685277/1125899906842624;
    dft(i) = 5034666847464581/140737488355328 - (4661696384331003*x(i))/8796093022208;
    dfs(i) = - (2412762607603083*x(i))/8796093022208 - 8094894005075091/140737488355328;
end

%%
% % plot TSR and SSR against FAL
figure;
hold on
t = plot(x, tsr, 'b.');
s = plot(x, ssr, 'r.');
plot(fit_x, fit_ft, 'b--')
plot(fit_x, fit_fs, 'r--')
hold off
ylabel('net solar radiation (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
title('net SW flux over fal at (75, 75), May-Jul','FontSize',24);
legend([t,s],{'TSR','SSR'},'FontSize',20)
set(gca,'FontSize',20)
set(gcf, 'Position', [10, 10, 650, 400]);

%%
% % plot dTSR/dFAL and dSSR/dFAL
figure
hold on
t=plot(x,dft,'b.--');
s=plot(x,dfs,'r.--');
hold off
ylabel('\partial F/\partial \alpha (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
title('Derivative of net SW flux over fal at (75, 75), May-Jul','FontSize',24);
legend([t,s],{'TSR','SSR'},'FontSize',20)
set(gca,'FontSize',20)
set(gcf, 'Position', [10, 10, 650, 400]);
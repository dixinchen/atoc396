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

% get albedo data
data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
x = zeros(4,38);
for i=1:38
    for j = 4:7
        x(j-3,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j));
    end
end 
x = x(:);


% get SSR data
data = ncread([gen_path,'data_1_m/GSSR_1_m38.nc'],'ssr');
ssr = zeros(4,38);
for i=1:38
    for j = 4:7
        ssr(j-3,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
    end
end
ssr = ssr(:);

% get TSR data
data = ncread([gen_path,'data_1_m/GTSR_1_m38.nc'],'tsr');
tsr = zeros(4,38);
for i=1:38
    for j = 4:7
        tsr(j-3,i)=squeeze(data(idx_lon_nn, idx_lat_nn, (i-1)*12+j))/(24*3600);
    end
end
tsr = tsr(:);

% % plot TSR and SSR against FAL
figure;
hold on
t = plot(x, tsr, '.');
s = plot(x, ssr, '.');
hold off
ylabel('net solar radiation (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
title('net SW flux over fal at (75, 75), Apr-Jul','FontSize',24);
legend([t,s],{'TSR','SSR'},'FontSize',20)
set(gca,'FontSize',20)
set(gcf, 'Position', [10, 10, 650, 400]);
saveas(gcf,strcat(figure_path, "amjj.png"));


% get climatological TISR data
data = ncread([gen_path,'data_1_m/insolation_clim_mean.nc'],'ALLSKY_TOA_SW_DWN');
data = squeeze(data.*(3600000/(24*3600)));
tisr = zeros(4,38);
for i = 1:38
    for j = 4:7
        tisr(j-3,i) = data(j);
    end
end
tisr = tisr(:);


ref = zeros(size(tsr));
tr = zeros(size(tsr));
for i=1:size(tsr)
    al = x(i);
    fs = ssr(i);
    ft = tsr(i);
    s = tisr(i);
    ref(i) = (- al^2*s^2 + ft*al^2*s + al*fs^2 + 2*al*s^2 - 2*ft*al*s - s^2 + ft*s)/(al^2*fs^2 - al^2*s^2 + 2*al*s^2 - s^2);
    tr(i) = -(fs*s + al*fs*ft - 2*al*fs*s - al^2*fs*ft + al^2*fs*s)/(al^2*fs^2 - al^2*s^2 + 2*al*s^2 - s^2);
end

dft_dal = zeros(size(tsr));
dfs_dal = zeros(size(tsr));
for i=1:size(tsr)
    al = x(i);
    s = tisr(i);
    r = ref(i);
    t = tr(i);
    dft_dal(i) = (s*t^2)/(al*r - 1) - (al*r*s*t^2)/(al*r - 1)^2;
    dfs_dal(i) = (s*t)/(al*r - 1) - (r*s*t*(al - 1))/(al*r - 1)^2;
end

% plot dft_dal and dfs_dal over fal
figure
hold on
t = plot(x, dft_dal, '.');
s = plot(x, dfs_dal, '.');
hold off
ylabel('\partial F / \partial\alpha (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
title('Change of net SFC flux w.r.t. albedo, Apr-Jul','FontSize',24);
legend([t,s],{'TOA','SFC'},'FontSize',20,'Location','southeast')
set(gca,'FontSize',20)
set(gcf, 'Position', [10, 10, 650, 400]);
saveas(gcf,strcat(figure_path, "diff_amjj.png"));
%%
data = ncread([gen_path,'data_1_m/Test2007_2016/dfal_1_m10.nc'],'d_fal');
x = squeeze(data(idx_lon_nn, idx_lat_nn, 4:7,:));
da = x(:);

dft = dft_dal((4*28)+1:end).* da;
dfs = dfs_dal((4*28)+1:end).* da;

fit_x = linspace(max(da), min(da), 1000);
fit_dft = polyval(polyfit(da, dft, 2),fit_x);
fit_dfs = polyval(polyfit(da, dfs, 2),fit_x);

figure
hold on
s = plot(da, dfs, 'r.');
t = plot(da, dft, 'b.');
plot(fit_x, fit_dft, 'b--')
plot(fit_x, fit_dfs, 'r--')
hold off
ylabel('\partial F (W/m^2)','FontSize',24);xlabel('dfal (%)','FontSize',24);
title('Theoretical \Deltaflux due to albedo, Apr-Jul','FontSize',24);
legend([s,t],{'SFC','TOA'},'FontSize',20)
set(gca,'FontSize',20)
% set(gcf, 'Position', [10, 10, 650, 400]);
saveas(gcf,strcat(figure_path, "the_fal_feedback_anomaly_amjj.png"));
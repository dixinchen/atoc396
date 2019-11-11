clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';

m = ["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"];

trend = "deseasonal";
% trend = "interannual";

% var = "wv";
var = "fal";

month = 6;
num = 120;

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
clim_mean = ncread([gen_path, 'data_1_m/Test2007_2016/fal_mean_10.nc'],'fal_mean');
clim_mean = squeeze(clim_mean(idx_lon_nn, idx_lat_nn, month));


x = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
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
set(gca,'FontSize',20)
ylabel('dRa (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
title(strcat(m(month)," 2007-2016, ", trend),'FontSize',24);
if trend == "deseasonal"
    legend([k,nn,mean],{'Kernel','NN','Clim. mean'},'FontSize',20);
end
if trend == "interannual"
    legend([k,nn],{'Kernel','NN'},'FontSize',20);
end
saveas(gcf,figure_path+trend+'-drx-dx-'+var+month+'.png');

%% plot dra/a 12*10 (all months in 10 years)

x = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
if trend == "deseasonal"
    x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+1:end));
end
if trend == "interannual"
    x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+1:end));
    x = x.';
    X = 1:num;
    p = polyfit(X',x',1);
    x_dt = x'-p(1)*X'-p(2);
    x = x_dt;
end

dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
if trend == "deseasonal"
    dr_k = squeeze(dr(2,1,idx_lon_k, idx_lat_k,:,:));
    dr_k = dr_k(:);
end
if trend == "interannual"
    dr = squeeze(dr(2,1,idx_lon_k, idx_lat_k,:,:));
    dr = dr.';
    dr = dr(:).';
    X = 1:num;
    p = polyfit(X',dr',1);
    dr_dt = dr'-p(1)*X'-p(2);
    dr_k = dr_dt;
end

dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
if trend == "deseasonal"
    dr_nn = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,:,:));
    dr_nn = dr_nn(:);
end
if trend == "interannual"
    dr = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,:,:));
    dr = dr.';
    dr = dr(:).';
    X = 1:num;
    p = polyfit(X',dr',1);
    dr_dt = dr'-p(1)*X'-p(2);
    dr_nn = dr_dt;
end

figure;
hold on
k = plot(x,dr_k, '.','MarkerSize',10);
nn = plot(x,dr_nn,'.','MarkerSize',10);
hold off
set(gca,'FontSize',20)
ylabel('dRa (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
title(strcat("10 yrs deseasonal dR_a / fal at (75)"),'FontSize',20);
legend([k,nn],{'Kernel','NN'},'FontSize',20);
% saveas(gcf,figure_path+'-drx-dx-'+var+'-all-mths.png');
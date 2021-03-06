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
x = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
if trend == "deseasonal"
%     x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
    x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+1:end));
%     x = squeeze(x(idx_lon_nn, idx_lat_nn, month:12:12*37+month));
end
if trend == "interannual"
    x = squeeze(x(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month));
%     x = squeeze(x(idx_lon_nn, idx_lat_nn, month:12:12*37+month));
    x = x.';
    X = 1:num;
    p = polyfit(X',x',1);
    x_dt = x'-p(1)*X'-p(2);
    x = x_dt;
end

% get SSR data
dr = ncread([gen_path,'data_1_m/GSSR_1_m38.nc'],'ssr');
if trend == "deseasonal"
%     ssr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month))./(24*3600);
    ssr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+1:end))./(24*3600);
%     ssr = squeeze(dr(idx_lon_nn, idx_lat_nn, month:12:12*37+month))./(24*3600);
%     ssr = squeeze(dr(idx_lon_nn, idx_lat_nn, :))./2.628e+6;
end
if trend == "interannual"
    dr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month))./(24*3600);
%     dr = squeeze(dr(idx_lon_nn, idx_lat_nn, month:12:12*37+month))./(24*3600);
%     dr = squeeze(dr(idx_lon_nn, idx_lat_nn, :))./2.628e+6;
    dr = dr.';
%     X = 1:456;
    X = 1:num;
    p = polyfit(X',dr',1);
    ssr = dr'-p(1)*X'-p(2);
end

% get TSR data
dr = ncread([gen_path,'data_1_m/GTSR_1_m38.nc'],'tsr');
if trend == "deseasonal"
%     tsr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month))./(24*3600);
    tsr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+1:end))./(24*3600);
%     tsr = squeeze(dr(idx_lon_nn, idx_lat_nn, month:12:12*37+month))./(24*3600);
%     tsr = squeeze(dr(idx_lon_nn, idx_lat_nn, :))./2.628e+6;
end
if trend == "interannual"
    dr = squeeze(dr(idx_lon_nn, idx_lat_nn, 12*28+month:12:12*37+month))./(24*3600);
%     dr = squeeze(dr(idx_lon_nn, idx_lat_nn, month:12:12*37+month))./(24*3600);
%     dr = squeeze(dr(idx_lon_nn, idx_lat_nn, :))./2.628e+6;
    dr = dr.';
%     X = 1:456;
    X = 1:num;
    p = polyfit(X',dr',1);
    tsr = dr'-p(1)*X'-p(2);
end

% get TISR data
% tisr = zeros(10,1);
% for i = 2007:2016
%     dr = ncread(strcat(gen_path,'data_1_m/tisr2007_2016/tisr',num2str(i),'.nc'),'tisr');
%     tisr(i-2006) = squeeze(dr(idx_lon_nn, idx_lat_nn, month))./(24*3600);
% end

tisr = zeros(12,10);
for i = 2007:2016
    dr = ncread(strcat(gen_path,'data_1_m/tisr2007_2016/tsi',num2str(i),'.nc'),'TSI');
%     tisr(:,i-2006) = squeeze(dr(idx_lon_nn, idx_lat_nn, :))./(24*3600);
    tisr(:,i-2006) = squeeze(dr(:));
end
tisr = tisr(:);


%%% plot TSR and SSR against FAL
% figure;
% hold on
% t = plot(x, tsr, '.');
% s = plot(x, ssr, '.');
% hold off
% ylabel('net solar radiation (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
% title('net SW flux over fal at (75, 75) in Jun*10','FontSize',24);
% legend([t,s],{'TSR','SSR'},'FontSize',20)
% set(gca,'FontSize',20)


%%% calculate A (absorptance of atmos) and R (reflectance at TOA)
% syms al fs ft s r t
% eqns = [fs==(1-al)*s*(t/(1-al*r)), ft==s-s*r-s*al*t^2/(1-al*r)];
% S = solve(eqns,[r t]);

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

%%% diff
syms al fs ft s r t
ft = s-s*r-s*al*t^2/(1-al*r);
fs = (1-al)*s*(t/(1-al*r));
diff(ft,al)
diff(fs,al)


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

% % get NN method dra data
% dr = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
% dr_k = squeeze(dr(2,1,idx_lon_k, idx_lat_k,month,:));
% dr_k = dr_k(:);
% 
% dr = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
% dr_nn = squeeze(dr(2,1,idx_lon_nn, idx_lat_nn,month,:));
% dr_nn = dr_nn(:);

% fit_x = linspace(max(x), min(x), 1000);
% fit_dr = polyval(polyfit(x, dfs_dal, 2),fit_x);

figure
hold on
me = semilogy(x, dfs_dal, 'b.');
% plot(fit_x,fit_dr,'b--');
% k = plot(x, dr_k, '.');
% nn = plot(x, dr_nn, '.');
hold off
ylabel('\partial F_{SFC} / \partial\alpha (W/m^2)','FontSize',24);xlabel('fal (%)','FontSize',24);
% title('Change of net SFC flux w.r.t. albedo in June','FontSize',24);
% legend([me,k,nn],{'The','Kernel','NN'},'FontSize',20)
set(gca,'FontSize',20)
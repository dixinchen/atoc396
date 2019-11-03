clc;clear;close all;tic;

% num = 12*20;
gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test2007_2016/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';

%% set month and coordinate of interest
month = 1;
lon = 205;
lat = 70;
W = cos(lat'*pi/180);

% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'dtcwv_1_m10.nc'],'latitude');
idx_lon_nn = find(lon_nn==lon);
idx_lat_nn = find(lat_nn==lat);

% coord index for 2.5*2.5 resolution
lon_k = ncread([result_path,'Feedback_kernel_cld.nc'],'longitude');
lat_k = ncread([result_path,'Feedback_kernel_cld.nc'],'latitude');
idx_lon_k = find(lon_k==lon);
idx_lat_k = find(lat_k==lat);

num = 10;
dR_nn_m = zeros(12,10);
dR_k_m = zeros(12,10);
dfal_m = zeros(12,10);

for i = 1:12
    data = ncread([gen_path,'dfal_1_m10.nc'],'d_fal');
    data = squeeze(data(idx_lon_nn, idx_lat_nn, i, :));
    data = data.*W;
    data = data.';
%     X = 1:num;
%     p = polyfit(X',data',1);
%     data_dt = data'-p(1)*X'-p(2);
%     data = data_dt.';
    data = 2*(data-repmat(min(data),[1 num]))./repmat(max(data)-min(data),[1 num])-1;
    dfal_m(i,:) = data;
    
    data = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
    data = squeeze(data(2,1,idx_lon_k, idx_lat_k, i, :));
    data = data.*W;
    data = data.';
%     X = 1:num;
%     p = polyfit(X',data',1);
%     data_dt = data'-p(1)*X'-p(2);
%     data = data_dt.';
    data = 2*(data-repmat(min(data),[1 num]))./repmat(max(data)-min(data),[1 num])-1;
    dR_k_m(i,:) = data;
    
    data = ncread([result_path,'Feedback_NN_cld.nc'],'dR');
    data = squeeze(data(2,1,idx_lon_nn, idx_lat_nn, i, :));
    data = data.*W;
    data = data.';    
%     X = 1:num;
%     p = polyfit(X',data',1);
%     data_dt = data'-p(1)*X'-p(2);
%     data = data_dt.';
    data = 2*(data-repmat(min(data),[1 num]))./repmat(max(data)-min(data),[1 num])-1;
    dR_nn_m(i,:) = data;
end

%% plot
% figure;
% hold on
% for i = (8:11)
%     plot(dfal_m(i,:),dR_nn_m(i,:),'x','MarkerSize',14, 'DisplayName',['NN-',num2str(i)]);
%     plot(dfal_m(i,:),dR_k_m(i,:),'.','MarkerSize',14,'DisplayName',['K-',num2str(i)]);   
% end
% legend show


%%
sel_month = [8 9 10 11];
% sel_month = [5 6];
% sel_month = [1 11 12];
% sel_month = [1:12];
dfal_sel_m = dfal_m(sel_month,:);
dfal_sel_m = dfal_sel_m(:);
dR_nn_sel_m = dR_nn_m(sel_month,:);
dR_nn_sel_m = dR_nn_sel_m(:);
dR_k_sel_m = dR_k_m(sel_month,:);
dR_k_sel_m = dR_k_sel_m(:);

% dfal_sel_m = dfal_sel_m(dfal_sel_m > -1e-3 & dfal_sel_m < 1e-3);
% dR_nn_sel_m = dR_nn_sel_m(dR_nn_sel_m < 0.3 & dR_nn_sel_m > -0.1);
% dR_k_sel_m = dR_k_sel_m(dR_k_sel_m < 0.3 & dR_k_sel_m > -0.1);


% fit
fit_data = linspace(max(dfal_sel_m), min(dfal_sel_m), 1000);
fit_dR_k_sel_m = polyval(polyfit(dfal_sel_m, dR_k_sel_m, 1),fit_data);
fit_dR_nn_sel_m = polyval(polyfit(dfal_sel_m, dR_nn_sel_m, 2),fit_data);

figure;
hold on
h1 = plot(dfal_sel_m,dR_k_sel_m,'.','MarkerSize',15);%, 'DisplayName','Kernel'
h2 = plot(dfal_sel_m,dR_nn_sel_m,'.','MarkerSize',15);%, 'DisplayName','NN'
h3 = plot(fit_data,fit_dR_k_sel_m,'b--','LineWidth',2);%'LineWidth',8,
h4 = plot(fit_data,fit_dR_nn_sel_m,'r--', 'LineWidth',2);%'LineWidth',8,
ylabel('dRa [W/m^2]','FontSize',16);xlabel('dfal','FontSize',16);title('1997-2016 JJA','FontSize',16)
legend([h1 h2],{'Kernel','NN'},'FontSize',16)

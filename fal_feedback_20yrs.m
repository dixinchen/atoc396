clc;clear;close all;tic;

% num = 12*20;
gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';

%% set month and coordinate of interest

lon = 205;
lat = 70;
W = cos(lat'*pi/180);

% coord index for 2.5*2.5 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test1997_2016/dtcwv_1_m20.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test1997_2016/dtcwv_1_m20.nc'],'latitude');
idx_lon_nn = find(lon_nn==lon);
idx_lat_nn = find(lat_nn==lat);

num = 20;
dR_m = zeros(12,20);
dfal_m = zeros(12,20);
fal_m = zeros(12,20);

for i = 1:12
    data = ncread([gen_path,'data_1_m/Test1997_2016/dfal_1_m20.nc'],'d_fal');
    dfal = squeeze(data(idx_lon_nn, idx_lat_nn, i, :));
    dfal = reshape(dfal,[1 1 num]).*repmat(W,[1 1 num]);
    dfal = squeeze(nansum(dfal,1));
    dfal = dfal.';
    X = 1:num;
    p = polyfit(X', dfal',1);
    dfal_dt = dfal'-p(1)*X'-p(2);
    dfal_dt = dfal_dt.';
%     dfal_dt = 2*(dfal_dt-repmat(min(dfal_dt),[1 num]))./repmat(max(dfal_dt)-min(dfal_dt),[1 num])-1;
    dfal_m(i,:) = dfal_dt;
    
    fal = ncread([gen_path, 'data_1_m/GFAL_1_m38.nc'], 'fal');
    fal = fal(idx_lon_nn, idx_lat_nn, 12*18+i:12:12*37+i);
    fal = reshape(fal,[1 1 20]).*repmat(W,[1 1 20]);
    fal = squeeze(nansum(fal,1));
    fal_m(i,:) = fal.';
    
    dR = ncread([result_path,'Feedback_NN_cld_20yrs.nc'],'dR');
    dR = squeeze(dR(1, 1, idx_lon_nn, idx_lat_nn, i, :));
    dR_nn = reshape(dR,[1 1 num]).*repmat(W,[1 1 num]);
    dR_nn = squeeze(nansum(dR_nn,1));
    dR_nn = dR_nn.';
%     X = 1:num;
%     p = polyfit(X',dR_nn',1);
%     dR_nn_dt = dR_nn'-p(1)*X'-p(2);
%     dR_nn_dt = dR_nn_dt.';
%     dR_nn_dt = 2*(dR_nn_dt-repmat(min(dR_nn_dt),[1 num]))./repmat(max(dR_nn_dt)-min(dR_nn_dt),[1 num])-1;
%     dR_m(i,:) = dR_nn_dt;

    dR_nn = 2*(dR_nn-repmat(min(dR_nn),[1 num]))./repmat(max(dR_nn)-min(dR_nn),[1 num])-1;   
    dR_m(i,:) = dR_nn;
end

month = 9;

clim_m = ncread([gen_path,'data_1_m/Test1997_2016/fal_mean_20.nc'],'fal_mean');
clim_m = squeeze(clim_m(idx_lon_nn, idx_lat_nn, month));

% % load EARi data (fal monthly mean)
% fal = ncread([gen_path, 'data_1_m/GFAL_1_m38.nc'], 'fal');
% fal = fal(idx_lon_nn, idx_lat_nn, 12*18+month:12:12*37+month);
% fal = reshape(fal,[1 1 20]).*repmat(W,[1 1 20]);
% fal = squeeze(nansum(fal,1));
% fal = fal.';


%% plot every month
figure;
hold on
for i = [1:12]
    a = plot(fal_m(i,:), dR_m(i,:),'.','MarkerSize',10,'DisplayName',['month: ',num2str(i)]);
end
xline(clim_m);
hold on
legend('FontSize',14,'Location','SouthEast')
ylabel('dRa [W/m^2]','FontSize',16);xlabel('fal','FontSize',16);title('1997-2016 NN estimation','FontSize',16)


%%
% select month
sel_month = [1:12];
dfal_sel_m = dfal_m(sel_month,:);
dfal_sel_m = dfal_sel_m(:);
dR_sel_m = dR_m(sel_month,:);
dR_sel_m = dR_sel_m(:);

%% zooming in
dfal_sel_m = dfal_sel_m(dfal_sel_m > -1e-3 & dfal_sel_m < 1e-3);
dR_nn_sel_m = dR_nn_sel_m(dR_nn_sel_m < 0.3 & dR_nn_sel_m > -0.1);
dR_k_sel_m = dR_k_sel_m(dR_k_sel_m < 0.3 & dR_k_sel_m > -0.1);

%%
% fit
fit_data = linspace(max(dfal_sel_m), min(dfal_sel_m), 1000);
fit_dR_nn_sel_m = polyval(polyfit(dfal_sel_m, dR_sel_m, 2),fit_data);

figure;
hold on
h2 = plot(dfal_sel_m,dR_sel_m,'.','MarkerSize',15);%, 'DisplayName','NN'
h4 = plot(fit_data,fit_dR_nn_sel_m,'r--', 'LineWidth',2);%'LineWidth',8,
ylabel('dRa [W/m^2]','FontSize',16);xlabel('dfal','FontSize',16);title('1997-2016 JJA','FontSize',16)
legend([h2],{'NN'},'FontSize',16)
%%
% figure;
% % plot(fal, dtsr, '--');
% % hold on
% % plot(fal, dR_k, '.');
% % hold on
% plot(dfal_dt, -dR_nn_dt, '.');
% hold on
% % plot(dfal, -dR_nn, '.');
% % hold on
% % xline(cli_mean(month));
% % hold on
% % legend('EARi', 'Kernel', 'NN', 'Clim. mean');

%% load climatological mean
% cli_mean = ncread([gen_path, 'data_1_m/Test1997_2016/fal_mean_20.nc'],'fal_mean');
% cli_mean = cli_mean(idx_lon_nn, idx_lat_nn, :).*W;

%% load EARi data (fal monthly mean)
% fal = ncread([gen_path, 'data_1_m/GFAL_1_m38.nc'], 'fal');
% fal = fal(idx_lon_nn, idx_lat_nn, 12*18+month:12:12*37+month);
% fal = reshape(fal,[1 1 20]).*repmat(W,[1 1 20]);
% fal = squeeze(nansum(fal,1));
% fal = fal.';

%% load true data
% 
% % TOA SW
% dR = ncread([gen_path,'data_1_m/Test1997_2016/dtsr_1_m20.nc'],'d_tsr');
% dR = squeeze(dR(idx_lon_nn, idx_lat_nn, month, :));
% dtsr = reshape(dR,[1 1 20]).*repmat(W,[1 1 20]);
% dtsr = squeeze(nansum(dtsr,1));
% dtsr = dtsr.';
% 
% % SFC SW
% dR = ncread([gen_path,'data_1_m/Test1997_2016/dssr_1_m20.nc'],'d_ssr');
% dR = squeeze(dR(idx_lon_nn, idx_lat_nn, month, :));
% dssr = reshape(dR,[1 1 20]).*repmat(W,[1 1 20]);
% dssr = squeeze(nansum(dssr,1));
% dssr = dssr.';
% X = 1:20;
% p = polyfit(X', dssr',1);
% dssr_dt = dssr'-p(1)*X'-p(2);
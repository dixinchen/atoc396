% (2000-1979)*12+4:12:(2015-1979)*12+4
clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
my_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-20/';

lat1 = 70;
lat_mid = 80;
lat2 = 90;
lon1 = [180 270 0 900];
% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'latitude');
idx_lon1 = zeros(1,4);
% for i =1:4
%     idx_lon1(i) = find(lon_nn==lon1(i));
% end
idx_lat1_nn = find(lat_nn==lat1);
idx_lat2_nn = find(lat_nn==lat2);
idx_lat_mid = find(lat_nn==lat_mid);

la1 = 69.5;
la2 = 89.5;
lon2 = [0.5, 90.5, 180.5, 270.5];
lon_ebaf = ncread([my_path,'swd_sfc.nc'],'lon');
lat_ebaf = ncread([my_path,'swd_sfc.nc'],'lat');
idx_lon2 = zeros(1,4);
for i =1:4
    idx_lon2(i) = find(lon_ebaf==lon2(i));
end
idx_lat1_ebaf = find(lat_ebaf==la1);
idx_lat2_ebaf = find(lat_ebaf==la2);

lon_tmp_nn = lon_nn(:);
lat_tmp_nn = lat_nn(idx_lat2_nn:idx_lat1_nn);
%% albedo anomaly (Delta alpha)
data = ncread([gen_path,'data_1_m/anomaly/dfal_1_m38.nc'],'d_fal');
dfal = data(:,idx_lat2_nn:idx_lat1_nn,4:9,2000-1979+1:2015-1979+1);

%% albedo monthly mean (alpha)
fal = NaN(360,21,6,16);
for month = 4:9
    data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
    fal(:,:,month-3,:) = data(:,idx_lat2_nn:idx_lat1_nn,(2000-1979)*12+month:12:(2015-1979)*12+month);
end
%% dASR_atm = [all sky - clear sky] net toa solar radiation anomaly
data = ncread([gen_path,'data_1_m/anomaly/dtsr_1_m38.nc'],'d_tsr');
datac = ncread([gen_path,'data_1_m/anomaly/dtsrc_1_m38.nc'],'d_tsrc');
dasr_atm = data(:,idx_lat2_nn:idx_lat1_nn,4:9,2000-1979+1:2015-1979+1)-datac(:,idx_lat2_nn:idx_lat1_nn,4:9,2000-1979+1:2015-1979+1);


%% SW downward radiation monthly mean (SWD_SFC)
% EBAF data
data = ncread([my_path,'swd_sfc.nc'],'sfc_sw_down_all_mon');
swd = NaN(360,21,6,16);
for i=1:16
    for j=4:9
        swd(:,:,j-3,i) = squeeze(data(:, 160:180, 12*(i-1)+(j-2)));
    end
end

% swd_sfc clim. mean
data = ncread([my_path,'swf_sfc_clim.nc'],'sfc_sw_down_all_clim');
swd_clim = NaN(360,21,6);
for j = 4:9
    swd_clim(:,:,j-3) = data(:,160:180, j);
end

% calculate swd_sfc anomaly dswd
dswd = NaN(360,21,6,16);
for i = 1:16
    for j = 1:6
        dswd(:,:,j,i)=swd(:,:,j,i)-swd_clim(:,:,j);
    end
end


%% calculate contributions
con_atm = dasr_atm+(1-fal).*dswd;
% con_sfc = (1-dfal).*swd;
con_sfc = -dfal.*swd;
con = con_atm+con_sfc;

con_atm(:,:,3,2) = NaN(size(con_atm,1),size(con_atm,2));
con_sfc(:,:,3,2) = NaN(size(con_sfc,1),size(con_sfc,2));

%% cos lat weighting

w = zeros(4,90,11);
wp = zeros(1,6);

lon = lon_nn(:);
lat = lat_nn(idx_lat2_nn:idx_lat1_nn);
[~,y]=meshgrid(lon,lat);
w_arctic = cos(y'*pi/180);
wp(1,1) = sum(w_arctic,[1 2]);

idx = [181 271 1 91 270 360 90 180];
% cos lat weighting
for i=1:4
    lon = lon_nn(idx(i):idx(i+4));
    lat = lat_nn(idx_lat_mid:idx_lat1_nn);
    [~,y]=meshgrid(lon,lat);
    w(i,:,:) = cos(y'*pi/180);
    wp(1,i+1) = sum(w(i,:,:),[2 3]);
end

lon = lon_nn(:);
lat = lat_nn(idx_lat2_nn:idx_lat_mid);
numlon = size(lon,1);
numlat = size(lat,1);
[~,y]=meshgrid(lon,lat);
w_region5 = cos(y'*pi/180);
wp(1,end) = sum(w_region5,[1 2]);


%% weighted mean in six regions

means = zeros(6,3,6,16); % Dimensions: '6'-regions; '3'-CON/CON_SFC/CON_ATM; '6'-months; '16'-yrs

% 1. arctic region
for month = 1:6
    for yr = 1:16
        % CON TOTAL
        data = squeeze(con(:,idx_lat2_nn:idx_lat1_nn,month,yr));
        data = data.*w_arctic;
        data = data(:);
        data = sum(data);
        data = data/wp(1);
        means(1,1,month,yr)=data;
        % CON SFC
        data = squeeze(con_sfc(:,idx_lat2_nn:idx_lat1_nn,month,yr));
        data = data.*w_arctic;
        data = data(:);
        data = sum(data);
        data = data/wp(1);
        means(1,2,month,yr)=data;
        % CON ATM
        data = squeeze(con_atm(:,idx_lat2_nn:idx_lat1_nn,month,yr));
        data = data.*w_arctic;
        data = data(:);
        data = sum(data);
        data = data/wp(1);
        means(1,3,month,yr)=data;
    end
end

% 2. region 1-4

for month = 1:6
    for i = 1:4 
        for yr = 1:16
            data = squeeze(con(idx(i):idx(i+4),idx_lat_mid:idx_lat1_nn,month,yr));
            data = data.*w(i);
            data = data(:);
            data = sum(data);
            data = data/wp(i+1);
            means(i+1,1,month,yr)=data;
        end
        for yr = 1:16
            data = squeeze(con_sfc(idx(i):idx(i+4),idx_lat_mid:idx_lat1_nn,month,yr));
            data = data.*w(i);
            data = data(:);
            data = sum(data);
            data = data/wp(i+1);
            means(i+1,2,month,yr)=data;
        end
        for yr = 1:16
            data = squeeze(con_atm(idx(i):idx(i+4),idx_lat_mid:idx_lat1_nn,month,yr));
            data = data.*w(i);
            data = data(:);
            data = sum(data);
            data = data/wp(i+1);
            means(i+1,3,month,yr)=data;
        end
    end
end

% 3. region 5
for month = 1:6
    for yr = 1:16
        % CON TOTAL
        data = squeeze(con(:,idx_lat2_nn:idx_lat_mid,month,yr));
        data = data.*w_region5;
        data = data(:);
        data = sum(data);
        data = data/wp(6);
        means(6,1,month,yr)=data;
        % CON SFC
        data = squeeze(con_sfc(:,idx_lat2_nn:idx_lat_mid,month,yr));
        data = data.*w_region5;
        data = data(:);
        data = sum(data);
        data = data/wp(6);
        means(6,2,month,yr)=data;
        % CON ATM
        data = squeeze(con_atm(:,idx_lat2_nn:idx_lat_mid,month,yr));
        data = data.*w_region5;
        data = data(:);
        data = sum(data);
        data = data/wp(6);
        means(6,3,month,yr)=data;
    end
end


%% calculate variance (along years)
% Dimensions: '6'-regions; '3'-CON/CON_SFC/CON_ATM; '6'-months

var = nanvar(means,0,4);

%% plot
axes = [0.5 6.5 0 100; 0.5 6.5 0 100; 0.5 6.5 0 250; 0.5 6.5 0 250; 0.5 6.5 0 100; 0.5 6.5 0 100];
m = ["Apr" "May" "Jun" "Jul" "Aug" "Sept"];
names = {'Arctic'; 'Region I'; 'Region II'; 'Region III'; 'Region IV'; 'Region V'};
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for month=1:6  
    subplot(3,2,month);
    b = bar(squeeze(var(:,:,month)));
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'b';
    b(3).FaceColor = 'r';
    axis(axes(month,:))
    set(gca,'xtick',1:6,'xticklabel',names)
    set(gca,'FontSize',18)
    if month ==1
        legend([b(1),b(2),b(3)],{'total ASR','sfc.','atm.'},'FontSize',18)
    end
end


%% save data (N/A)
% filename = [my_path,'var_zhan.nc'];
% if exist(filename,'file')
%     delete(filename);
% end
% %save dim
% struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
% netcdf_write(filename,lon_nn,struct_tmp);
% struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
% netcdf_write(filename,lat_nn(idx_lat2_nn:idx_lat1_nn),struct_tmp);
% %save data
% struct_tmp = struct('name','var_total','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Total variance',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status1] = netcdf_write(filename,var_total,struct_tmp);
% 
% struct_tmp = struct('name','var_atm','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','ATM contribution to dASR',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status2] = netcdf_write(filename,var_atm,struct_tmp);
% 
% struct_tmp = struct('name','var_sfc','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','SFC contribution to dASR',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status3] = netcdf_write(filename,var_sfc,struct_tmp);
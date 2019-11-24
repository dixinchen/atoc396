% (2000-1979)*12+4:12:(2015-1979)*12+4
clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
my_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-20/';

lat1 = 30;
lat2 = 90;
% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'latitude');
idx_lat1_nn = find(lat_nn==lat1);
idx_lat2_nn = find(lat_nn==lat2);

la1 = 29.5;
la2 = 89.5;
lon_ebaf = ncread([my_path,'swd_sfc.nc'],'lon');
lat_ebaf = ncread([my_path,'swd_sfc.nc'],'lat');
idx_lat1_ebaf = find(lat_ebaf==la1);
idx_lat2_ebaf = find(lat_ebaf==la2);

%% albedo anomaly (Delta alpha)
data = ncread([gen_path,'data_1_m/anomaly/dfal_1_m38.nc'],'d_fal');
dfal = data(:,idx_lat2_nn:idx_lat1_nn,4:9,2000-1979+1:2015-1979+1);

%% albedo monthly mean (alpha)
fal = NaN(360,idx_lat1_nn-idx_lat2_nn+1,6,16);
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
swd = NaN(360,idx_lat1_nn-idx_lat2_nn+1,6,16);
for i=1:16
    for j=4:9
        swd(:,:,j-3,i) = squeeze(data(:, idx_lat2_ebaf:-1:idx_lat1_ebaf, 12*(i-1)+(j-2)));
    end
end

% swd_sfc clim. mean
data = ncread([my_path,'swf_sfc_clim.nc'],'sfc_sw_down_all_clim');
swd_clim = NaN(360,idx_lat1_nn-idx_lat2_nn+1,6);
for j = 4:9
    swd_clim(:,:,j-3) = data(:, idx_lat2_ebaf:-1:idx_lat1_ebaf, j);
end

% calculate swd_sfc anomaly dswd
dswd = NaN(360,idx_lat1_nn-idx_lat2_nn+1,6,16);
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

%% save data before computing var of regional mean
% calculate var along years in each grid
gridded_var = zeros(3,6,360,idx_lat1_nn-idx_lat2_nn+1);
for month = 1:6
    for i=1:360
        for j=1:idx_lat1_nn-idx_lat2_nn+1
            current_con = squeeze(con(i,j,month,:));
            gridded_var(1,month,i,j)=nanvar(current_con);
            current_con = squeeze(con_sfc(i,j,month,:));
            gridded_var(2,month,i,j)=nanvar(current_con);
            current_con = squeeze(con_atm(i,j,month,:));
            gridded_var(3,month,i,j)=nanvar(current_con);        
        end
    end
end

% save data (June only)
filename = [my_path,'var_zhan.nc'];
if exist(filename,'file')
    delete(filename);
end
%save dim
struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lon_nn,struct_tmp);
struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lat_nn(idx_lat2_nn:idx_lat1_nn),struct_tmp);
%save data
struct_tmp = struct('name','var_total','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Total variance',...
    'units','(W/m^2)^2','missing_values',-9999);
[status1] = netcdf_write(filename,gridded_var(1,3,:,:),struct_tmp);

struct_tmp = struct('name','var_sfc','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','ATM contribution to dASR',...
    'units','(W/m^2)^2','missing_values',-9999);
[status2] = netcdf_write(filename,gridded_var(2,3,:,:),struct_tmp);

struct_tmp = struct('name','var_atm','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','SFC contribution to dASR',...
    'units','(W/m^2)^2','missing_values',-9999);
[status3] = netcdf_write(filename,gridded_var(3,3,:,:),struct_tmp);
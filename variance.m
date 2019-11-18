clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';
my_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/';

m = ["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"];

trend = "deseasonal";
% trend = "interannual";
var = "fal";
month = 6;
num = 10;
% lon = 75; lat = 75;
lat1 = 70;
lat2 = 90;

% coord index for 2.5*2.5 resolution
lon_k = ncread([result_path,'Feedback_kernel_cld.nc'],'longitude');
lat_k = ncread([result_path,'Feedback_kernel_cld.nc'],'latitude');
idx_lat1_k = find(lat_k==lat1);
idx_lat2_k = find(lat_k==lat2);

% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'latitude');
idx_lat1_nn = find(lat_nn==lat1);
idx_lat2_nn = find(lat_nn==lat2);

lon_tmp_nn = lon_nn(:);
lat_tmp_nn = lat_nn(idx_lat2_nn:idx_lat1_nn);
lon_tmp_k = lon_k(:);
lat_tmp_k = lat_k(idx_lat2_k:idx_lat1_k);

% get NN dra data (TOA SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRa');
dra_nn = squeeze(data(1,1,:,idx_lat2_nn:idx_lat1_nn, month, :));
% dra_nn = squeeze(data(1,1,:,:, month, :));
var_dra_nn = nanvar(dra_nn,0,3);

% get NN drc data (TOA SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRc');
drc_nn = squeeze(data(1,1,:,idx_lat2_nn:idx_lat1_nn, month, :));
var_drc_nn = nanvar(drc_nn,0,3);

% get kernel dra data (TOA SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRa');
dra_k = squeeze(data(1,1,:,idx_lat2_k:idx_lat1_k, month, :));
var_dra_k = nanvar(dra_k,0,3);

% get kernel drc data (TOA SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRc');
drc_k = squeeze(data(1,1,:,idx_lat2_k:idx_lat1_k, month, :));
var_drc_k = nanvar(drc_k,0,3);

% get NN dra data (SFC SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRa');
dra_nn = squeeze(data(2,1,:,idx_lat2_nn:idx_lat1_nn, month, :));
svar_dra_nn = nanvar(dra_nn,0,3);

% get NN drc data (SFC SW)
data = ncread([result_path, 'Feedback_NN_cld.nc'],'dRc');
drc_nn = squeeze(data(2,1,:,idx_lat2_nn:idx_lat1_nn, month, :));
svar_drc_nn = nanvar(drc_nn,0,3);

% get kernel dra data (SFC SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRa');
dra_k = squeeze(data(2,1,:,idx_lat2_k:idx_lat1_k, month, :));
svar_dra_k = nanvar(dra_k,0,3);

% get kernel drc data (SFC SW)
data = ncread([result_path, 'Feedback_kernel_cld.nc'],'dRc');
drc_k = squeeze(data(2,1,:,idx_lat2_k:idx_lat1_k, month, :));
svar_drc_k = nanvar(drc_k,0,3);


% filename = [my_path,'var_dra_drc_globle_nn.nc'];
% if exist(filename,'file')
%     delete(filename);
% end
% %save dim
% struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
% % netcdf_write(filename,lon_tmp_nn,struct_tmp);
% netcdf_write(filename,lon_nn,struct_tmp);
% struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
% netcdf_write(filename,lat_nn,struct_tmp);
% %save data
% struct_tmp = struct('name','var_toa_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA NN dRa',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status1] = netcdf_write(filename,var_dra_nn,struct_tmp);
% 
% struct_tmp = struct('name','var_toa_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA NN dRc',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status2] = netcdf_write(filename,var_drc_nn,struct_tmp);
% 
% struct_tmp = struct('name','var_sfc_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA NN dRa',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status3] = netcdf_write(filename,svar_dra_nn,struct_tmp);
% 
% struct_tmp = struct('name','var_sfc_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA NN dRc',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status4] = netcdf_write(filename,svar_drc_nn,struct_tmp);
% 
% 
% 
% filename = [my_path,'var_dra_drc_globle_kernel.nc'];
% if exist(filename,'file')
%     delete(filename);
% end
% %save dim
% struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
% netcdf_write(filename,lon_k,struct_tmp);
% struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
% netcdf_write(filename,lat_k,struct_tmp);
% %save data
% struct_tmp = struct('name','var_toa_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA kernel dRa',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status1] = netcdf_write(filename,var_dra_k,struct_tmp);
% 
% struct_tmp = struct('name','var_toa_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA kernel dRc',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status2] = netcdf_write(filename,var_drc_k,struct_tmp);
% 
% struct_tmp = struct('name','var_sfc_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA kernel dRa',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status3] = netcdf_write(filename,svar_dra_k,struct_tmp);
% 
% struct_tmp = struct('name','var_sfc_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
% struct_tmp.dim = {'longitude';'latitude'};
% struct_tmp.att = struct('short_name','Variance TOA kernel dRc',...
%     'units','(W/m^2)^2','missing_values',-9999);
% [status4] = netcdf_write(filename,svar_drc_k,struct_tmp);
% 
% toc;

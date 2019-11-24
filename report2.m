clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
my_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-26/';

% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'latitude');
% 
% data = ncread([gen_path, 'data_1_m/Test2007_2016/dtsr_1_m10.nc'], "d_ssr");
% x = squeeze(data(:,:,6,:));
% interan_x = zeros(size(x));
% for i=1:size(x,1)
%     for j = 1:size(x,2)
%         a = squeeze(x(i,j,:));
%         a = a';
%         X = 1:10;
%         p=polyfit(X',a',1);
%         interan_x(i,j,:)=a'-p(1)*X'-p(2);
%     end
% end
% var = nanvar(interan_x,0,3);
% var_with_lt = nanvar(x,0,3);

cloud = ["hcc" "mcc" "lcc"];
% cloud = ["tciw" "tclw"];
x = zeros(size(cloud,1),size(lon_nn,1),size(lat_nn,1),10);
for i=size(cloud,1)
    data = ncread(strcat(gen_path,'data_1_m/Test2007_2016/d',cloud(i),'_1_m10.nc'),strcat("d_",cloud(i)));
    x(i,:,:,:) = squeeze(data(:,:,6,:));
end
y = squeeze(nansum(x,1));
var_with_lt = nanvar(y,0,3);

%%
filename = [my_path,'var_cloud_deseasonal_nov23.nc'];
if exist(filename,'file')
    delete(filename);
end
struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lon_nn,struct_tmp);
struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lat_nn,struct_tmp);
struct_tmp = struct('name','var_tsr','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance of TSR',...
    'units','(W/m^2)^2','missing_values',-9999);
[status1] = netcdf_write(filename,var_with_lt,struct_tmp);
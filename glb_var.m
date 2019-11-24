clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';
my_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-20/';

var = "fal";
month = 6;
num = 10;

% coord index for 2.5*2.5 resolution
lon_k = ncread([result_path,'Feedback_kernel_cld.nc'],'longitude');
lat_k = ncread([result_path,'Feedback_kernel_cld.nc'],'latitude');

% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'latitude');
%%
var_nn = zeros(2,3,size(lon_nn,1),size(lat_nn,1));
name = ["dRa","dRc","dRw"];
for toa_sfc = 1:2
    for variables = 1:3
        data = ncread([result_path, 'Feedback_NN_cld.nc'],name(variables));
        x = squeeze(data(toa_sfc,1,:,:, month, :));
        ita_x = zeros(size(x));
        for i=1:size(x,1)
            for j = 1:size(x,2)
                a = squeeze(x(i,j,:));
                a = a';
                X = 1:10;
                p=polyfit(X',a',1);
                ita_x(i,j,:)=a'-p(1)*X'-p(2);
            end
        end
        var_nn(toa_sfc,variables,:,:) = nanvar(ita_x,0,3);
    end
end
toc;

%%
var_k = zeros(2,3,size(lon_k,1),size(lat_k,1));
name = ["dRa","dRc","dRw"];
for toa_sfc = 1:2
    for variables = 1:3
        data = ncread([result_path, 'Feedback_kernel_cld.nc'],name(variables));
        x = squeeze(data(toa_sfc,1,:,:, month, :));
        ita_x = zeros(size(x));
        for i=1:size(x,1)
            for j = 1:size(x,2)
                a = squeeze(x(i,j,:));
                a = a';
                X = 1:10;
                p=polyfit(X',a',1);
                ita_x(i,j,:)=a'-p(1)*X'-p(2);
            end
        end
        var_k(toa_sfc,variables,:,:) = nanvar(ita_x,0,3);
    end
end
toc;

%save NN data
filename = [my_path,'var_global_june_nn.nc'];
if exist(filename,'file')
    delete(filename);
end

%save dim
struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lon_nn,struct_tmp);
struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lat_nn,struct_tmp);
%save data
struct_tmp = struct('name','var_toa_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance TOA NN dRa',...
    'units','(W/m^2)^2','missing_values',-9999);
[status1] = netcdf_write(filename,squeeze(var_nn(1,1,:,:)),struct_tmp);

struct_tmp = struct('name','var_toa_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance TOA NN dRc',...
    'units','(W/m^2)^2','missing_values',-9999);
[status2] = netcdf_write(filename,squeeze(var_nn(1,2,:,:)),struct_tmp);

struct_tmp = struct('name','var_toa_drw','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance TOA NN dRw',...
    'units','(W/m^2)^2','missing_values',-9999);
[status3] = netcdf_write(filename,squeeze(var_nn(1,3,:,:)),struct_tmp);

struct_tmp = struct('name','var_sfc_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance SFC NN dRa',...
    'units','(W/m^2)^2','missing_values',-9999);
[status4] = netcdf_write(filename,squeeze(var_nn(2,1,:,:)),struct_tmp);

struct_tmp = struct('name','var_sfc_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance SFC NN dRc',...
    'units','(W/m^2)^2','missing_values',-9999);
[status5] = netcdf_write(filename,squeeze(var_nn(2,2,:,:)),struct_tmp);

struct_tmp = struct('name','var_sfc_drw','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance SFC NN dRw',...
    'units','(W/m^2)^2','missing_values',-9999);
[status6] = netcdf_write(filename,squeeze(var_nn(2,3,:,:)),struct_tmp);


% save kernel data
filename = [my_path,'var_global_june_kernel.nc'];
if exist(filename,'file')
    delete(filename);
end
%save dim
struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lon_k,struct_tmp);
struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lat_k,struct_tmp);
%save data
struct_tmp = struct('name','var_toa_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance TOA kernel dRa',...
    'units','(W/m^2)^2','missing_values',-9999);
[status1] = netcdf_write(filename,squeeze(var_k(1,1,:,:)),struct_tmp);


struct_tmp = struct('name','var_toa_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance TOA kernel dRc',...
    'units','(W/m^2)^2','missing_values',-9999);
[status2] = netcdf_write(filename,squeeze(var_k(1,2,:,:)),struct_tmp);

struct_tmp = struct('name','var_toa_drw','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance TOA kernel dRw',...
    'units','(W/m^2)^2','missing_values',-9999);
[status3] = netcdf_write(filename,squeeze(var_k(1,3,:,:)),struct_tmp);

struct_tmp = struct('name','var_sfc_dra','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance SFC kernel dRa',...
    'units','(W/m^2)^2','missing_values',-9999);
[status4] = netcdf_write(filename,squeeze(var_k(2,1,:,:)),struct_tmp);

struct_tmp = struct('name','var_sfc_drc','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance SFC kernel dRc',...
    'units','(W/m^2)^2','missing_values',-9999);
[status5] = netcdf_write(filename,squeeze(var_k(2,2,:,:)),struct_tmp);

struct_tmp = struct('name','var_sfc_drw','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'longitude';'latitude'};
struct_tmp.att = struct('short_name','Variance SFC kernel dRw',...
    'units','(W/m^2)^2','missing_values',-9999);
[status6] = netcdf_write(filename,squeeze(var_k(2,3,:,:)),struct_tmp);
%calculate the [20]-year mean
clc;clear;close all;tic;
rad = {'ssrc';'ssr';'tsrc';'tsr';'strc';'str';'ttrc';'ttr'};
data_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/G';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test1997_2016/';
% vars = {'tcc';'tcwv';'sp';'tco3';'fal';'tciw';'tclw';'hcc';'mcc';'lcc';'skt';'t500';'t200';'t10';'q700';'q500';'q200'};
vars = {'tcwv';'sp';'tco3';'fal';'tciw';'tclw';'hcc';'mcc';'lcc';'skt';'t500';'t200';'t10';'q700';'q500';'q200'};
nvar = size(vars,1);
unit_v = {'kg/m2';'Pa';'kg/m2';'0,1';'kg/m2';'kg/m2';'';'';'';'K';'K';'K';'K';'kg/kg';'kg/kg';'kg/kg'};
lon = ncread([data_path,'SSRC_1_m38.nc'],'longitude');
lat = ncread([data_path,'SSRC_1_m38.nc'],'latitude');
numlon = size(lon,1);
numlat = size(lat,1);
addpath('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test1997_2016/');

numy = 20;

for i = 1:8
    filename = [result_path,rad{i},'_mean_',num2str(numy),'.nc'];
    if exist(filename,'file')
    delete(filename);
    end
    month = 1:12;
    struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
    netcdf_write(filename,lon,struct_tmp);
    struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
    netcdf_write(filename,lat,struct_tmp);
    struct_tmp = struct('name','month','type','dim','nc_type','NC_SHORT');
    netcdf_write(filename,month,struct_tmp);
     X = ncread([data_path,upper(rad{i}),'_1_m38.nc'],rad{i});
     X = X(:,:,12*(38-numy)+1:end);
     X = X/24/3600;
    X = reshape(X,numlon,numlat,12,numy);
    Rmean = squeeze(nanmean(X,4));
    
    struct_tmp = struct('name',[rad{i},'_mean'],'type','var','nc_type','NC_FLOAT','var_name',[]);
    struct_tmp.dim = {'longitude';'latitude';'month'};
    struct_tmp.att = struct('short_name',['The monthly mean of ',upper(rad{i}),' for the period from ',num2str(2016-numy+1),' to 2016'],...
        'units',unit_v{i},'missing_values',-9999);
    [status1] = netcdf_write(filename,Rmean,struct_tmp);
    toc;
end

for i = 1:nvar
    filename = [result_path, vars{i},'_mean_',num2str(numy),'.nc'];
    if exist(filename,'file')
    delete(filename);
    end
    month = 1:12;
    struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
    netcdf_write(filename,lon,struct_tmp);
    struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
    netcdf_write(filename,lat,struct_tmp);
    struct_tmp = struct('name','month','type','dim','nc_type','NC_SHORT');
    netcdf_write(filename,month,struct_tmp);
     X = ncread([data_path,upper(vars{i}),'_1_m38.nc'],vars{i});
     X = X(:,:,12*(38-numy)+1:end);
     
    X = reshape(X,numlon,numlat,12,numy);
    Rmean = squeeze(nanmean(X,4));
    
    struct_tmp = struct('name',[vars{i},'_mean'],'type','var','nc_type','NC_FLOAT','var_name',[]);
    struct_tmp.dim = {'longitude';'latitude';'month'};
    struct_tmp.att = struct('short_name',['The monthly mean of ',upper(vars{i}),' for the period from ',num2str(2016-numy+1),' to 2016'],...
        'units',unit_v{i},'missing_values',-9999);
    [status1] = netcdf_write(filename,Rmean,struct_tmp);
    toc;
end
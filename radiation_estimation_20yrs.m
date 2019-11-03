%testing the performance of NN model on radiation flux estimation
clc;clear; close all ;tic;
data_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/';
model_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/Model_simplied/';
var_position = {'toa';'sfc'};
var_sky = {'all';'clr'};
var_wv ={'SW';'LW'};
%% determine the input variables
lon = ncread([data_path,'GSSR_1_m38.nc'],'longitude');
lat = ncread([data_path,'GSSR_1_m38.nc'],'latitude');

numlon = size(lon,1);
numlat = size(lat,1);
numy = 20;

Lon = zeros(numlon,2);
Lon(:,1) = [0:180 179:-1:1]';
Lon(:,2) = sin(lon*pi/180);
[x,y]=meshgrid(Lon(:,1),lat);%add longitude and latitude as inputs
W = cos(y'*pi/180);
Wper = nansum(nansum(W));
Lonc = zeros(numlon*numlat,2);
Lonc(:,1) = reshape(x',numlon*numlat,1);
[x,y]=meshgrid(Lon(:,2),lat);%add longitude and latitude as inputs
Lonc(:,2) = reshape(x',numlon*numlat,1);
Lonc1 = repmat(Lonc,[numy,1]);
Lat = reshape(y',numlon*numlat,1);
Lat1 = repmat(Lat,[numy,1]);

num = numlon*numlat*numy;
RF = nan(2,2,2,numlon,numlat,12,numy);
BR = nan(2,2,2,numlon,numlat,12,numy);
for var_1 = 1:2
    var1 = var_position{var_1};
    for var_2 = 1:2
        var2 = var_sky{var_2};
        for var_3 = 1:2
            var3 = var_wv{var_3};
            if strcmp(var1,'toa')
                if strcmp(var2,'all')
                    if strcmp(var3,'SW')
                        rad = 'tsr';
                        vars = {'tciw';'tclw';'hcc';'mcc';'lcc';'tcwv';'fal';'sp';'tco3'};
                        sw_flag = 1;
                    elseif strcmp(var3,'LW')
                        rad = 'ttr';
                        vars = {'skt';'t500';'t200';'t10';'q700';'q500';'q200';'hcc';'mcc';'lcc'};
                        sw_flag = 0;
                    end
                elseif strcmp(var2,'clr')
                    if strcmp(var3,'SW')
                        rad = 'tsrc';
                        vars = {'tcwv';'sp';'tco3';'fal'};
                        sw_flag = 1;
                    elseif strcmp(var3,'LW')
                        rad = 'ttrc';
                        vars = {'skt';'t500';'t200';'t10';'q700';'q500';'q200'};
                        sw_flag = 0;
                    end
                end
            elseif strcmp(var1,'sfc')
                if strcmp(var2,'all')
                    if strcmp(var3,'SW')
                        rad = 'ssr';
                        vars = {'tciw';'tclw';'hcc';'mcc';'lcc';'tcwv';'fal';'sp';'tco3'};
                        sw_flag = 1;
                    elseif strcmp(var3,'LW')
                        rad = 'str';
                        vars = {'skt';'t500';'t200';'t10';'tcwv';'hcc';'mcc';'lcc'};
                        sw_flag = 0;
                    end
                elseif strcmp(var2,'clr')
                    if strcmp(var3,'SW')
                        rad = 'ssrc';
                        vars = {'tcwv';'sp';'tco3';'fal'};
                        sw_flag = 1;
                    elseif strcmp(var3,'LW')
                        rad = 'strc';
                        vars = {'skt';'t500';'t200';'t10';'tcwv'};
                        sw_flag = 0;
                    end
                end
            end
            
            
            nvar = size(vars,1);
            %Calculate the last ten years from 2007 to 2016
            Rad = ncread([data_path,'G',upper(rad),'_1_m38.nc'],rad);
            Rad = Rad/24/3600;
            Rad = Rad(:,:,18*12+1:end);
            Rad = reshape(Rad,[numlon,numlat,12,20]);
            
            Rad_nn = zeros(numlon,numlat,12,20);
            load([model_path,'Net_',upper(rad),'.mat']);
            for i = 1:12
                SamIn = zeros(num,nvar+3);
                for k = 1:nvar
                    X = ncread([data_path,'G',upper(vars{k}),'_1_m38.nc'],vars{k});
                    X = X(:,:,18*12+1:end);
                    X = X(:,:,i:12:end);
                    SamIn(:,k) = reshape(X,[num,1]);
                end
                SamIn(:,nvar+1:nvar+3) = [Lonc1 Lat1];
                OutNN = NN_model(SamIn,Net,i);
%                 load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
%                 OutNN = sim(net,double(SamIn'))';
                if sw_flag
                    OutNN = OutNN.*(OutNN>0);
                end
                Rad_nn(:,:,i,:) = reshape(OutNN,[numlon,numlat,1,20]);
            end
            RF(var_1,var_2,var_3,:,:,:,:) = Rad_nn;
            BR(var_1,var_2,var_3,:,:,:,:) = Rad_nn-Rad;
            toc;
        end
    end
end
addpath('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/');
month = 1:12;
year = 1997:2016;

filename = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/Radiation_estimation_20yrs.nc';
if exist(filename,'file')
    delete(filename);
end

%save dim
struct_tmp = struct('name','TOA_SFC','type','dim','nc_type','NC_BYTE');
netcdf_write(filename,[1 2],struct_tmp);
struct_tmp = struct('name','All_Clr','type','dim','nc_type','NC_BYTE');
netcdf_write(filename,[1 2],struct_tmp);
struct_tmp = struct('name','SW_LW','type','dim','nc_type','NC_BYTE');
netcdf_write(filename,[1 2],struct_tmp);

struct_tmp = struct('name','longitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lon,struct_tmp);
struct_tmp = struct('name','latitude','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,lat,struct_tmp);
struct_tmp = struct('name','month','type','dim','nc_type','NC_FLOAT');
netcdf_write(filename,month,struct_tmp);
struct_tmp = struct('name','year','type','dim','nc_type','NC_SHORT');
netcdf_write(filename,year,struct_tmp);
%save data
struct_tmp = struct('name','R','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'All_Clr';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Radiation Flux calculated by NN Model',...
    'units','W/m2','missing_values',-9999);
[status1] = netcdf_write(filename,RF,struct_tmp);

struct_tmp = struct('name','B','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'All_Clr';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Radiation bias (Rnn-ERAi)',...
    'units','W/m2','missing_values',-9999);
[status2] = netcdf_write(filename,BR,struct_tmp);
toc;

function [ status ] = netcdf_write(filename,data,structure,if_disp)
if ~exist('if_disp','var') || isempty(if_disp)
    if_disp = 0;
end

if exist(filename,'file')
    if if_disp, disp(['Appending to - ',filename]); end
    nc = netcdf.open(filename,'NC_WRITE');
    netcdf.reDef(nc);
else
    if if_disp, disp(['Creating - ',filename]); end
    nc = netcdf.create(filename,'NC_WRITE');
end

if ~exist('structure','var')
    structure.name = inputname(2);
    structure.type = 'var';
end

if ~isfield(structure,'name')
    structure.name = inputname(2);
end
if ~isfield(structure,'type')
    structure.type = 'var';
end
if ~isfield(structure,'nc_type')
    if isinteger(data)
        structure.nc_type = 'NC_INT';
    elseif isfloat(data)
        structure.nc_type = 'NC_DOUBLE';
    elseif isnumeric(data)
        structure.nc_type = 'NC_DOUBLE';
    elseif isstr(data)
        structure.nc_type = 'NC_CHAR';
    elseif islogical(data)
        structure.nc_type = 'NC_BYTE';
    end
end
if strcmpi(structure.nc_type,'NC_BYTE')
    data = int8(data);
elseif strcmpi(structure.nc_type,'NC_CHAR')
    data = char(data);
elseif strcmpi(structure.nc_type,'NC_SHORT')
    data = int16(data);
elseif strcmpi(structure.nc_type,'NC_INT')
    data = int32(data);
elseif strcmpi(structure.nc_type,'NC_FLOAT')
    data = single(data);
elseif strcmpi(structure.nc_type,'NC_DOUBLE')
    data = double(data);
else
    disp(['ERROR: unrecognizable nc_type!']);
    netcdf.close(nc);
    status = -1;
    return;
end

if strcmpi(structure.type,'variable') ||  strcmpi(structure.type,'var')
    
    % dimension
    if isfield(structure,'dim')
        nc_dim = getfield(structure,'dim');
    elseif isfield(structure,'dimension')
        nc_dim = getfield(structure,'dimension');
    else
        for idim = 1:ndims(data)
            nc_dim(idim) = size(data,idim);
        end
    end
    %
    if ~iscell(nc_dim) && isstr(nc_dim)
        nc_dim = cellstr(nc_dim);
    end
    %
    if isnumeric(nc_dim) % new dimensions
        if prod(nc_dim) ~= length(data(:))
            disp('ERROR: assigned dimensions NOT matching data!');
            netcdf.close(nc);
            status = -1;
            return;
        end
        % def. dimensions
        for idim = 1:length(nc_dim)
            dim_data(idim) = netcdf.defDim(nc,[structure.name,'_dim',num2str(idim)],nc_dim(idim));
        end
    else % existent dimensions
        total_size = 1;
        for idim = 1:length(nc_dim)
            dim_data(idim) = netcdf.inqDimID(nc,nc_dim{idim});
        end
    end
    % def. variable
    var_data = netcdf.defVar(nc,structure.name,structure.nc_type,[dim_data]);
    % put variable attributes
    if isfield(structure,'attribute')
        nc_att = getfield(structure,'attribute');
    elseif isfield(structure,'attr')
        nc_att = getfield(structure,'attr');
    elseif isfield(structure,'att')
        nc_att = getfield(structure,'att');
    else
        nc_att = [];
    end
    if ~isempty(nc_att)
        att_names = fieldnames(nc_att);
        for iatt = 1:length(att_names)
            netcdf.putAtt(nc,var_data,att_names{iatt},eval(['nc_att.',att_names{iatt}]));
        end
        if isfield(nc_att,'missing_value')
            data(isnan(data)) = nc_att.missing_value;
        end
    end
    %
    netcdf.endDef(nc);
    % put var
    netcdf.putVar(nc,var_data,data);
    
elseif strcmpi(structure.type,'dimension') ||  strcmpi(structure.type,'dim')
    
    % dimension
    if isfield(structure,'if_unlim') && structure.if_unlim == 1
        dim_data = netcdf.defDim(nc,structure.name,netcdf.getConstant('NC_UNLIMITED'));
    else
        dim_data = netcdf.defDim(nc,structure.name,length(data(:)));
    end
    % variable
    var_data = netcdf.defVar(nc,structure.name,structure.nc_type,[dim_data]);
    % put variable attributes
    if isfield(structure,'attribute')
        nc_att = getfield(structure,'attribute');
    elseif isfield(structure,'attr')
        nc_att = getfield(structure,'attr');
    elseif isfield(structure,'att')
        nc_att = getfield(structure,'att');
    else
        nc_att = [];
    end
    if ~isempty(nc_att)
        att_names = fieldnames(nc_att);
        for iatt = 1:length(att_names)
            netcdf.putAtt(nc,var_data,att_names{iatt},eval(['nc_att.',att_names{iatt}]));
        end
        if isfield(nc_att,'missing_value')
            data(isnan(data)) = nc_att.missing_value;
        end
    end
    %
    netcdf.endDef(nc);
    % put var
    netcdf.putVar(nc,var_data,data);
    
elseif strcmpi(structure.type,'attribute') ||  strcmpi(structure.type,'att')
    
    if isfield(structure,'var_name') % variable attribute
        var_data = netcdf.inqVarID(nc,var_name);
        netcdf.putAtt(nc,var_data,structure.name,data);
    else % global attribute
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),structure.name,data);
    end
    netcdf.endDef(nc);
    
else
    
    disp('ERROR: Unrecognizable type!');
    status = -1;
    
end

netcdf.close(nc);
status = 0;
end
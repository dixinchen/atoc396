%use NN model to calculate radiative feedback in all sky
clc;clear; close all ;tic;
data_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/';
model_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/Model_simplied/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
var_position = {'toa';'sfc'};
var_wv ={'SW';'LW'};
%% determine the input variables
lon = ncread([data_path,'GSSR_1_m38.nc'],'longitude');
lat = ncread([data_path,'GSSR_1_m38.nc'],'latitude');

numlon = size(lon,1);
numlat = size(lat,1);
ini_y = 1997;
fin_y = 2016;
numy = fin_y-ini_y+1;

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
RF_t = zeros(2,2,numlon,numlat,12,numy);
RF_w = zeros(2,2,numlon,numlat,12,numy);
RF_a = zeros(2,2,numlon,numlat,12,numy);
RF_c = zeros(2,2,numlon,numlat,12,numy);
RF_nn = zeros(2,2,numlon,numlat,12,numy);
RF_ano = zeros(2,2,numlon,numlat,12,numy);


for var_1 = 1:2
    var1 = var_position{var_1};
    for var_3 = 1:2
        var3 = var_wv{var_3};
        if strcmp(var1,'toa')
            if strcmp(var3,'SW')
               rad = 'tsr';
               vars = {'tciw';'tclw';'hcc';'mcc';'lcc';'tcwv';'fal';'sp';'tco3'};
                sw_flag = 1;
            elseif strcmp(var3,'LW')
                rad = 'ttr';
                vars = {'skt';'t500';'t200';'t10';'q700';'q500';'q200';'hcc';'mcc';'lcc'};
                sw_flag = 0;
            end
            
        elseif strcmp(var1,'sfc')
            
            if strcmp(var3,'SW')
                rad = 'ssr';
                vars = {'tciw';'tclw';'hcc';'mcc';'lcc';'tcwv';'fal';'sp';'tco3'};
                sw_flag = 1;
            elseif strcmp(var3,'LW')
                rad = 'str';
                vars = {'skt';'t500';'t200';'t10';'tcwv';'hcc';'mcc';'lcc'};
                sw_flag = 0;
            end
            
        end
        
        
        nvar = size(vars,1);
        %Calculate feedback from 1997 to 2016
        Radclim = zeros(numlon,numlat,12);
        load([model_path,'Net_',upper(rad),'.mat']);
        for i = 1:12
            SamIn = zeros(numlon*numlat,nvar+3);
            for k = 1:nvar
                X = ncread([data_path,'Test1997_2016/',vars{k},'_mean_20.nc'],[vars{k},'_mean']);
                SamIn(:,k) = reshape(X(:,:,i),[numlon*numlat,1]);
            end
            SamIn(:,nvar+1:nvar+3) = [Lonc Lat];
            OutNN = NN_model(SamIn,Net,i);
%             load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
%             OutNN = sim(net,double(SamIn'))';
            Radclim(:,:,i) = reshape(OutNN,[numlon,numlat,1]);
        end
        if sw_flag
            Radclim = Radclim.*(Radclim>0);
        end
        
        
        if strcmp(var3,'SW')
            %calculate water vapor and albedo
            Xx = ncread([data_path,'G',upper(vars{6}),'_1_m38.nc'],vars{6});
            Xx = Xx(:,:,(ini_y-1979)*12+1:(fin_y-1979+1)*12);
            dRx = zeros(numlon,numlat,12,numy);
            for i = 1:12
                SamIn = zeros(numlon*numlat,nvar+3);
                for k = 1:nvar
                    X = ncread([data_path,'Test1997_2016/',vars{k},'_mean_20.nc'],[vars{k},'_mean']);
                    SamIn(:,k) = reshape(X(:,:,i),[numlon*numlat,1]);
                end
                SamIn(:,nvar+1:nvar+3) = [Lonc Lat];
%                 load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
                for iy = 0:numy-1
                    X = Xx(:,:,iy*12+i);
                    SamIn(:,6) = reshape(X,[numlon*numlat,1]);
                    OutNN = NN_model(SamIn,Net,i);
%                     OutNN = sim(net,double(SamIn'))';
                    dRx(:,:,i,iy+1) = reshape(OutNN,[numlon,numlat,1,1]);
                end
            end
            dRx = dRx.*(dRx>0);
            RF_w(var_1,1,:,:,:,:) = dRx-repmat(Radclim,[1,1,1,numy]);
            %albedo
            Xx = ncread([data_path,'G',upper(vars{7}),'_1_m38.nc'],vars{7});
            Xx = Xx(:,:,(ini_y-1979)*12+1:(fin_y-1979+1)*12);
            dRx = zeros(numlon,numlat,12,numy);
            for i = 1:12
                SamIn = zeros(numlon*numlat,nvar+3);
                for k = 1:nvar
                    X = ncread([data_path,'Test1997_2016/',vars{k},'_mean_20.nc'],[vars{k},'_mean']);
                    SamIn(:,k) = reshape(X(:,:,i),[numlon*numlat,1]);
                end
                SamIn(:,nvar+1:nvar+3) = [Lonc Lat];
%                 load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
                for iy = 0:numy-1
                    X = Xx(:,:,iy*12+i);
                    SamIn(:,7) = reshape(X,[numlon*numlat,1]);
                    OutNN = NN_model(SamIn,Net,i);
%                     OutNN = sim(net,double(SamIn'))';
                    dRx(:,:,i,iy+1) = reshape(OutNN,[numlon,numlat,1,1]);
                end
            end
            dRx = dRx.*(dRx>0);
            RF_a(var_1,1,:,:,:,:) = dRx-repmat(Radclim,[1,1,1,numy]);
            
            
            %cloud feedback
            dRx = zeros(numlon,numlat,12,numy);
            for i = 1:12
                SamIn = zeros(num,nvar+3);
                for k = 1:nvar
                    if k<6
                        Xx = ncread([data_path,'G',upper(vars{k}),'_1_m38.nc'],vars{k});
                        Xx = Xx(:,:,(ini_y-1979)*12+1:(fin_y-1979+1)*12);
                        X = Xx(:,:,i:12:end);
%                         X = X(:,:,1:numy); 
                        SamIn(:,k) = reshape(X,[num,1]);
                    else
                        Xm = ncread([data_path,'Test1997_2016/',vars{k},'_mean_20.nc'],[vars{k},'_mean']);
                        X = reshape(Xm(:,:,i),[numlon*numlat,1]);
                        SamIn(:,k) = repmat(X,[numy,1]);
                    end
                end
                SamIn(:,nvar+1:nvar+3) = [Lonc1 Lat1];
                OutNN = NN_model(SamIn,Net,i);
%                 load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
%                 OutNN = sim(net,double(SamIn'))';
                dRx(:,:,i,:) = reshape(OutNN,[numlon,numlat,1,numy]);
            end
            dRx = dRx.*(dRx>0);
            RF_c(var_1,1,:,:,:,:) = dRx-repmat(Radclim,[1,1,1,numy]);
        elseif strcmp(var3,'LW')

            %calculate temperature
            dRx = zeros(numlon,numlat,12,numy);
            for i = 1:12
                SamIn = zeros(num,nvar+3);
                for k = 1:nvar
                    if k<5
                        Xx = ncread([data_path,'G',upper(vars{k}),'_1_m38.nc'],vars{k});
                        Xx = Xx(:,:,(ini_y-1979)*12+1:(fin_y-1979+1)*12);
                        X = Xx(:,:,i:12:end);
%                         X = X(:,:,1:numy);
                        SamIn(:,k) = reshape(X,[num,1]);
                    else
                        Xm = ncread([data_path,'Test1997_2016/',vars{k},'_mean_20.nc'],[vars{k},'_mean']);
                        X = reshape(Xm(:,:,i),[numlon*numlat,1]);
                        SamIn(:,k) = repmat(X,[numy,1]);
                    end
                end
                SamIn(:,nvar+1:nvar+3) = [Lonc1 Lat1];
                OutNN = NN_model(SamIn,Net,i);
%                 load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
%                 OutNN = sim(net,double(SamIn'))';
                dRx(:,:,i,:) = reshape(OutNN,[numlon,numlat,1,numy]);
            end
            RF_t(var_1,2,:,:,:,:) = dRx-repmat(Radclim,[1,1,1,numy]);

            %calculate water vapor
            dRx = zeros(numlon,numlat,12,numy);
            for i = 1:12
                SamIn = zeros(num,nvar+3);
                for k = 1:nvar
                    if (k>4)&&(k<(nvar-2))
                        Xx = ncread([data_path,'G',upper(vars{k}),'_1_m38.nc'],vars{k});
                        Xx = Xx(:,:,(ini_y-1979)*12+1:(fin_y-1979+1)*12);
                        X = Xx(:,:,i:12:end);
                        SamIn(:,k) = reshape(X,[num,1]);
                    else
                        Xm = ncread([data_path,'Test1997_2016/',vars{k},'_mean_20.nc'],[vars{k},'_mean']);
                        X = reshape(Xm(:,:,i),[numlon*numlat,1]);
                        SamIn(:,k) = repmat(X,[numy,1]);
                    end
                end
                SamIn(:,nvar+1:nvar+3) = [Lonc1 Lat1];
                OutNN = NN_model(SamIn,Net,i);
%                 load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
%                 OutNN = sim(net,double(SamIn'))';
                dRx(:,:,i,:) = reshape(OutNN,[numlon,numlat,1,numy]);
            end
            RF_w(var_1,2,:,:,:,:) = dRx-repmat(Radclim,[1,1,1,numy]);
            %cloud feedback
            dRx = zeros(numlon,numlat,12,numy);
            for i = 1:12
                SamIn = zeros(num,nvar+3);
                for k = 1:nvar
                    if k>(nvar-3)
                        Xx = ncread([data_path,'G',upper(vars{k}),'_1_m38.nc'],vars{k});
                        Xx = Xx(:,:,(ini_y-1979)*12+1:(fin_y-1979+1)*12);
                        X = Xx(:,:,i:12:end);
                        SamIn(:,k) = reshape(X,[num,1]);
                    else
                        Xm = ncread([data_path,'Test1997_2016/',vars{k},'_mean_20.nc'],[vars{k},'_mean']);
                        X = reshape(Xm(:,:,i),[numlon*numlat,1]);
                        SamIn(:,k) = repmat(X,[numy,1]);
                    end
                end
                SamIn(:,nvar+1:nvar+3) = [Lonc1 Lat1];
                OutNN = NN_model(SamIn,Net,i);
%                 load([model_path,upper(rad),'/Gnet_',num2str(i),'_',rad,'.mat']);
%                 OutNN = sim(net,double(SamIn'))';
                dRx(:,:,i,:) = reshape(OutNN,[numlon,numlat,1,numy]);
            end
            RF_c(var_1,2,:,:,:,:) = dRx-repmat(Radclim,[1,1,1,numy]);
        end
        
        toc;
        dR = ncread([data_path,'Test1997_2016/d',rad,'_1_m20.nc'],['d_',rad]);
        RF_ano(var_1,var_3,:,:,:,:) = dR(:,:,:,1:numy);
        R = ncread([result_path,'Radiation_estimation_20yrs.nc'],'R');
        R = squeeze(R(var_1,1,var_3,:,:,:,1:numy));
        RF_nn(var_1,var_3,:,:,:,:) = R-repmat(Radclim,[1,1,1,numy]);
    end
end

RF_sum = RF_t+RF_w+RF_a+RF_c;

addpath('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/');
month = 1:12;
% year = 1997:2016;
year = ini_y:(ini_y+numy-1);

filename = [result_path,'Feedback_NN_cld_20yrs.nc'];
if exist(filename,'file')
    delete(filename);
end
%save dim
struct_tmp = struct('name','TOA_SFC','type','dim','nc_type','NC_BYTE');
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
struct_tmp = struct('name','dR','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Radiation anomaly',...
    'units','W/m2','missing_values',-9999);
[status1] = netcdf_write(filename,RF_ano,struct_tmp);

struct_tmp = struct('name','dRnn','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Total climate feedback by NN model changing together',...
    'units','W/m2','missing_values',-9999);
[status2] = netcdf_write(filename,RF_nn,struct_tmp);

struct_tmp = struct('name','dRt','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Temperature feedback',...
    'units','W/m2','missing_values',-9999);
[status3] = netcdf_write(filename,RF_t,struct_tmp);

struct_tmp = struct('name','dRw','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Water vapor feedback',...
    'units','W/m2','missing_values',-9999);
[status4] = netcdf_write(filename,RF_w,struct_tmp);

struct_tmp = struct('name','dRa','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Surface albedo feedback',...
    'units','W/m2','missing_values',-9999);
[status5] = netcdf_write(filename,RF_a,struct_tmp);

struct_tmp = struct('name','dRc','type','var','nc_type','NC_FLOAT','var_name',[]);
struct_tmp.dim = {'TOA_SFC';'SW_LW';'longitude';'latitude';'month';'year'};
struct_tmp.att = struct('short_name','Cloud feedback',...
    'units','W/m2','missing_values',-9999);
[status6] = netcdf_write(filename,RF_c,struct_tmp);
toc;
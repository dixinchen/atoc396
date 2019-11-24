% (2000-1979)*12+4:12:(2015-1979)*12+4


clc;clear;close all;tic;

gen_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/';
result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
result_path20 = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report5/output/';
figure_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-12/';
my_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-20/';

num = 10;

% 2000 - 2015 interannual


lat1 = 70;
lat2 = 90;
lon1 = [0 90 180 270];
lon2 = [0.5, 90.5, 180.5, 270.5];
% coord index for 1*1 resolution
lon_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'longitude');
lat_nn = ncread([gen_path, 'data_1_m/Test2007_2016/dtcwv_1_m10.nc'],'latitude');
idx_lon1 = zeros(1,4);
for i =1:4
    idx_lon1(i) = find(lon_nn==lon1(i));
end
idx_lat1_nn = find(lat_nn==lat1);
idx_lat2_nn = find(lat_nn==lat2);

%%
num = 16;

lon = ncread('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test2007_2016/dskt_1_m10.nc','longitude');
lat = ncread('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test2007_2016/dskt_1_m10.nc','latitude');
numlon = size(lon,1);
numlat = size(lat,1);
[~,y]=meshgrid(lon,lat);
Wnn = cos(y'*pi/180);

%%
la1 = 69.5;
la2 = 89.5;
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

%% test section

data = ncread([gen_path,'data_1_m/anomaly/dfal_1_m38.nc'],'d_fal');
data = data(:,:,4:9,2000-1979+1:2015-1979+1);

dfal = zeros(360,21,6,16);
for month = 6
    for i=360
        for j=21
            a = squeeze(data(i,j,month,:));
            x = 1:16;
            p = polyfit(x',a,1);
            aa = a-p(1)*x'-p(2);
            for k=16
                dfal(i,j,month,k)=aa(k);
            end
        end
    end
end
%%
% GFAL_1_m38.nc
data = ncread([gen_path,'data_1_m/GFAL_1_m38.nc'],'fal');
fal = zeros(360,21,6,16);
for i=1:16
    for j = 4:9
        if i==2 && j==6
            fal(:,:,j-3,i)=NaN;
        end
        fal(:,:,j-3,i)=squeeze(data(:, idx_lat2_nn:idx_lat1_nn, (20+i)*12+j));
    end
end

%%
data = ncread([gen_path,'data_1_m/anomaly/dtsr_1_m38.nc'],'d_tsr');
datac = ncread([gen_path,'data_1_m/anomaly/dtsrc_1_m38.nc'],'d_tsrc');
dasr_atm = zeros(360,21,6,16);
for i=1:16
    for j = 4:9
        if i==2 && j==6
            dasr_atm(:,:,j-3,i)=NaN;
        end
        dasr_atm(:,:,j-3,i)=squeeze(data(:, idx_lat2_nn:idx_lat1_nn, (20+i)*12+j)-datac(:, idx_lat2_nn:idx_lat1_nn, (20+i)*12+j));
    end
end

%%
data = ncread([my_path,'swd_sfc.nc'],'sfc_sw_down_all_mon');
swd = zeros(360,21,6,16);
for i=1:16
    for j = 4:9
        if i==2 && j==6
            swd(:,:,j-3,i)=NaN;
        end
        swd(:,:,j-3,i)=squeeze(data(:, 160:180, 12*(i-1)+j));
    end
end


% swf_sfc_clim.nc
data = ncread([my_path,'swf_sfc_clim.nc'],'sfc_sw_down_all_clim');
swd_clim = zeros(360,21,6);
for j = 4:9
    swd_clim(:,:,j)=data(:, 160:180, j);
end

dswd = zeros(360,21,6,16);
for i = 1:16
    for j = 4:9
        dswd(:,:,j-3,i)=swd(:,:,j-3,i)-swd_clim(:,:,j);
    end
end

con_atm = zeros(360,21,6,16);
for i = 1:16
    for j = 4:9
        con_atm(:,:,j-3,i)=dasr_atm(:,:,j-3,i)+(1-fal(:,:,j-3,i)).*dswd(:,:,j-3,i);
    end
end

con_sfc = zeros(360,21,6,16);
for i = 1:16
    for j = 4:9
        con_atm(:,:,j-3,i)=dfal(:,:,j-3,i).*swd(:,:,j-3,i);
    end
end

con = con_atm+con_sfc;

var_con = nanvar(squeeze(con(:,:,1,:)),0,3);
var_con_atm = nanvar(squeeze(con_atm(:,:,1,:)),0,3);
var_con_sfc = nanvar(squeeze(con_sfc(:,:,1,:)),0,3);

%%
tiledlayout(3,2)
for month = 1:6
    ax = nexttile;
    arctic = []
    bar(ax,a)
end
% Top bar graph
ax1 = nexttile;
bar(ax1,y)

% Bottom bar graph
ax2 = nexttile;
bar(ax2,y,'stacked')
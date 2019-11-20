%Fig4
clc;clear;close all;tic;

num = 12*10;

lon = ncread('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test2007_2016/dskt_1_m10.nc','longitude');
lat = ncread('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test2007_2016/dskt_1_m10.nc','latitude');
numlon = size(lon,1);
numlat = size(lat,1);
[~,y]=meshgrid(lon,lat);
Wnn = cos(y'*pi/180);

dskt = ncread('/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test2007_2016/dskt_1_m10.nc','d_skt');
Wper = sum(sum(Wnn));
dT = reshape(dskt,[numlon numlat 120]).*repmat(Wnn,[1 1 120]);
dT = squeeze(nansum(dT,1));
dT = squeeze(nansum(dT,1))/Wper;
X = 1:120;
p = polyfit(X',dT',1);
dT_dt = dT'-p(1)*X'-p(2);

Gr = cell(4,1);
for i = 1:4
    Gr{i,1}(1:120,1) = 0;
end
rad = {'ssrc';'strc';'tsrc';'ttrc';'ssr';'str';'tsr';'ttr'};
for i = 1:8    
    drad = ncread(['/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/data_1_m/Test2007_2016/d',rad{i},'_1_m10.nc'],['d_',rad{i}]);
    drad = reshape(drad,[numlon numlat 120]).*repmat(Wnn,[1 1 120]);
    drad = squeeze(nansum(drad,1));
    drad = squeeze(nansum(drad,1))/Wper;
    p = polyfit(X',drad',1);
    drad_dt = drad'-p(1)*X'-p(2);
    Gr{floor((i-1)/2)+1,1}(1:120,1) = Gr{floor((i-1)/2)+1,1}(1:120,1)+drad_dt;
end


result_path = '/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/sept/report2/output_ztt/';
Gk = cell(4,1);  %save data of time series
Gnn = cell(4,1);
%get kernel results
lon_k = ncread([result_path,'Feedback_kernel_cld.nc'],'longitude');
lat_k = ncread([result_path,'Feedback_kernel_cld.nc'],'latitude');
[~,y]=meshgrid(lon_k,lat_k);
Wk = cos(y'*pi/180);

%%
dRt_k = ncread([result_path,'Feedback_kernel_clr.nc'],'dRt');
[~,Gtk_dt] = Map2TS(Wk,dRt_k);
Gk{1,1}(1:120,2) = Gtk_dt(:,2);
Gk{2,1}(1:120,2) = Gtk_dt(:,1);
dRw_k = ncread([result_path,'Feedback_kernel_clr.nc'],'dRw');
[~,Gwk_dt] = Map2TS(Wk,dRw_k);
Gk{1,1}(1:120,3) = Gwk_dt(:,2);
Gk{2,1}(1:120,3) = Gwk_dt(:,1);
dRa_k = ncread([result_path,'Feedback_kernel_clr.nc'],'dRa');
[~,Gak_dt] = Map2TS(Wk,dRa_k);
Gk{1,1}(1:120,4) = Gak_dt(:,2);
Gk{2,1}(1:120,4) = Gak_dt(:,1);
sumdR_k = dRt_k+dRw_k+dRa_k;
[~,Gk_dt] = Map2TS(Wk,sumdR_k);
Gk{1,1}(1:120,6) = Gk_dt(:,2);
Gk{2,1}(1:120,6) = Gk_dt(:,1);


dRt_k = ncread([result_path,'Feedback_kernel_cld.nc'],'dRt');
[~,Gtk_dt] = Map2TS(Wk,dRt_k);
Gk{3,1}(1:120,2) = Gtk_dt(:,2);
Gk{4,1}(1:120,2) = Gtk_dt(:,1);
dRw_k = ncread([result_path,'Feedback_kernel_cld.nc'],'dRw');
[~,Gwk_dt] = Map2TS(Wk,dRw_k);
Gk{3,1}(1:120,3) = Gwk_dt(:,2);
Gk{4,1}(1:120,3) = Gwk_dt(:,1);
dRa_k = ncread([result_path,'Feedback_kernel_cld.nc'],'dRa');
[~,Gak_dt] = Map2TS(Wk,dRa_k);
Gk{3,1}(1:120,4) = Gak_dt(:,2);
Gk{4,1}(1:120,4) = Gak_dt(:,1);
dRc_k = ncread([result_path,'Feedback_kernel_cld.nc'],'dRc');
[~,Gck_dt] = Map2TS(Wk,dRc_k);
Gk{3,1}(1:120,5) = Gck_dt(:,2);
Gk{4,1}(1:120,5) = Gck_dt(:,1);
sumdR_k = dRt_k+dRw_k+dRa_k+dRc_k;
[~,Gk_dt] = Map2TS(Wk,sumdR_k);
Gk{3,1}(1:120,6) = Gk_dt(:,2);
Gk{4,1}(1:120,6) = Gk_dt(:,1);


%%

dRt_k = ncread([result_path,'Feedback_kernel_clr.nc'],'dRt');
a = Map2TS(Wk,dRt_k);
a_toa = a(:,1);



%%
%get NN results
dRt_nn = ncread([result_path,'Feedback_NN_clr.nc'],'dRt');
[~,Gtnn_dt] = Map2TS(Wnn,dRt_nn);
Gnn{1,1}(1:120,2) = Gtnn_dt(:,2);
Gnn{2,1}(1:120,2) = Gtnn_dt(:,1);
dRw_nn = ncread([result_path,'Feedback_NN_clr.nc'],'dRw');
[~,Gwnn_dt] = Map2TS(Wnn,dRw_nn);
Gnn{1,1}(1:120,3) = Gwnn_dt(:,2);
Gnn{2,1}(1:120,3) = Gwnn_dt(:,1);
dRa_nn = ncread([result_path,'Feedback_NN_clr.nc'],'dRa');
[~,Gann_dt] = Map2TS(Wnn,dRa_nn);
Gnn{1,1}(1:120,4) = Gann_dt(:,2);
Gnn{2,1}(1:120,4) = Gann_dt(:,1);
sumdR_nn = dRt_nn+dRw_nn+dRa_nn;
[~,Gn_dt] = Map2TS(Wnn,sumdR_nn);
Gnn{1,1}(1:120,6) = Gn_dt(:,2);
Gnn{2,1}(1:120,6) = Gn_dt(:,1);

dRt_nn = ncread([result_path,'Feedback_NN_cld.nc'],'dRt');
[~,Gtnn_dt] = Map2TS(Wnn,dRt_nn);
Gnn{3,1}(1:120,2) = Gtnn_dt(:,2);
Gnn{4,1}(1:120,2) = Gtnn_dt(:,1);
dRw_nn = ncread([result_path,'Feedback_NN_cld.nc'],'dRw');
[~,Gwnn_dt] = Map2TS(Wnn,dRw_nn);
Gnn{3,1}(1:120,3) = Gwnn_dt(:,2);
Gnn{4,1}(1:120,3) = Gwnn_dt(:,1);
dRa_nn = ncread([result_path,'Feedback_NN_cld.nc'],'dRa');
[~,Gann_dt] = Map2TS(Wnn,dRa_nn);
Gnn{3,1}(1:120,4) = Gann_dt(:,2);
Gnn{4,1}(1:120,4) = Gann_dt(:,1);
dRc_nn = ncread([result_path,'Feedback_NN_cld.nc'],'dRc');
[~,Gcnn_dt] = Map2TS(Wnn,dRc_nn);
Gnn{3,1}(1:120,5) = Gcnn_dt(:,2);
Gnn{4,1}(1:120,5) = Gcnn_dt(:,1);
sumdR_nn = dRt_nn+dRw_nn+dRa_nn+dRc_nn;
[~,Gn_dt] = Map2TS(Wnn,sumdR_nn);
Gnn{3,1}(1:120,6) = Gn_dt(:,2);
Gnn{4,1}(1:120,6) = Gn_dt(:,1);

dR_nn = ncread([result_path,'Feedback_NN_clr.nc'],'dRnn');
[~,Gnn_dt] = Map2TS(Wnn,dR_nn);
Gnn{1,1}(1:120,8) = Gnn_dt(:,2);
Gnn{2,1}(1:120,8) = Gnn_dt(:,1);
dR_nn = ncread([result_path,'Feedback_NN_cld.nc'],'dRnn');
[~,Gnn_dt] = Map2TS(Wnn,dR_nn);
Gnn{3,1}(1:120,8) = Gnn_dt(:,2);
Gnn{4,1}(1:120,8) = Gnn_dt(:,1);
%%
Ddate = zeros(num,7);
Ddate = char(Ddate);
for i = 2007:2016
    for j = 1:12
        if j <10
            Ddate((i-2007)*12+j,:) = [num2str(i),'-0',num2str(j)];
        else
            Ddate((i-2007)*12+j,:) = [num2str(i),'-',num2str(j)];
        end
    end
end
DD = datenum(Ddate);
%%
dRt_k = Gk{4,1}(:,2);
dRwv_k = Gk{4,1}(:,3);
dRfal_k = Gk{4,1}(:,4);
dRcc_k = Gk{4,1}(:,5);
sumdR_k = Gk{4,1}(:,6);

dRt_nn = Gnn{4,1}(:,2);
dRwv_nn = Gnn{4,1}(:,3);
dRfal_nn = Gnn{4,1}(:,4);
dRcc_nn = Gnn{4,1}(:,5);
sumdR_nn = Gnn{4,1}(:,6);
dRa = Gr{4,1}(:,1); 
%%
figure;
subplot(2,1,1);
plot(DD,dRt_k,'b.-');
hold on
plot(DD,dRt_nn,'r.-');
hold on
plot(DD,zeros(num,1),'k-'); title('(a) Temperature','FontSize',16);
 format bank;
 haxis = gca;
 Dat = {'2007';'2008';'2009';'2010';'2011';'2012';'2013';'2014';'2015';'2016'};
set(haxis,'XTick',DD(1:12:end));
set(gca, 'XTicklabel',Dat);
set(haxis,'fontsize',16);
 axis([min(DD) max(DD) -3.5 3.5]);dateaxis('x');
 ylabel({'W/m^{2}'},'FontSize',16);xlabel('Date','FontSize',16);
 legend('Kernel','NN','Orientation','horizontal');%legend('boxoff');

subplot(2,1,2);
plot(DD,dRwv_k,'b.-');
hold on
plot(DD,dRwv_nn,'r.-');
 hold on
 plot(DD,zeros(num,1),'k-'); title('(b) Water vapor','FontSize',16);
 format bank;
 haxis = gca;
 set(haxis,'XTick',DD(1:12:end));
set(gca, 'XTicklabel',Dat);
set(haxis,'fontsize',16);
 axis([min(DD) max(DD) -2.5 2.5]);dateaxis('x');%
 ylabel({'W/m^{2}'},'FontSize',16);xlabel('Date','FontSize',16);
 legend('Kernel','NN','Orientation','horizontal');%legend('boxoff');
figure;
subplot(2,1,1);
plot(DD,dRfal_k,'b.-');
hold on
plot(DD,dRfal_nn,'r.-');
 hold on
 plot(DD,zeros(num,1),'k-'); title('(c) Surface albedo','FontSize',16);
 format bank;
 haxis = gca;
 set(haxis,'XTick',DD(1:12:end));
set(gca, 'XTicklabel',Dat);
set(haxis,'fontsize',16);
 axis([min(DD) max(DD) -2.5 2.5]);dateaxis('x');%
 ylabel({'W/m^{2}'},'FontSize',16);xlabel('Date','FontSize',16);
 legend('Kernel','NN','Orientation','horizontal');%legend('boxoff');

subplot(2,1,2);
plot(DD,dRcc_k,'b.-');
hold on
plot(DD,dRcc_nn,'r.-');
hold on
 plot(DD,zeros(num,1),'k-'); title('(d) Cloud','FontSize',16);
 format bank;
 haxis = gca;
 set(haxis,'XTick',DD(1:12:end));
set(gca, 'XTicklabel',Dat);
set(haxis,'fontsize',16);
 axis([min(DD) max(DD) -2.5 2.5]);dateaxis('x');%
 ylabel({'W/m^{2}'},'FontSize',16);xlabel('Date','FontSize',16);
 legend('Kernel','NN','Orientation','horizontal');%legend('boxoff');
 figure;
subplot(2,1,1);
plot(DD,sumdR_k,'b.-');
hold on
plot(DD,sumdR_nn,'r.-');
hold on
plot(DD,dRa,'g.-');
hold on
 plot(DD,zeros(num,1),'k-'); title('(e) Total radiation anomalies','FontSize',16);
 format bank;
 haxis = gca;
 set(haxis,'XTick',DD(1:12:end));
set(gca, 'XTicklabel',Dat);
set(haxis,'fontsize',16);
 axis([min(DD) max(DD) -3 3]);dateaxis('x');%
 ylabel({'W/m^{2}'},'FontSize',16);xlabel('Date','FontSize',16);
 legend('Kernel','NN','ERAi','Orientation','horizontal');%legend('boxoff');
toc;
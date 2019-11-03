%
function [Gx,Gx_dt] = Map2TS(W,dRx)         
Wper = sum(sum(W));
dRx = squeeze(nansum(dRx,2));	% LW+SW
dRx = permute(dRx,[2 3 4 5 1]);
dRx = dRx.*repmat(W,[1 1 12 10 2]);
dRx = squeeze(nansum(dRx,1));
Gx = squeeze(nansum(dRx,1))/Wper;
Gx = reshape(Gx,[12*10,2]); 
%detrend
X = 1:120;
p = polyfit(X',Gx(:,1),1);
Gx_dt(:,1) = Gx(:,1)-p(1)*X'-p(2);
p = polyfit(X',Gx(:,2),1);
Gx_dt(:,2) = Gx(:,2)-p(1)*X'-p(2);
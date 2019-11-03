%NN model application
%NN model parameters are in Net
%Net, cell, month*parameter
%parameter: input range, Wij, bj, Wj, c, target range 

function SamOut = NN_model(SamIn,Net,i)
num = size(SamIn,1);
%input normalization
MinSam = Net{i,1}(:,1);
MaxSam = Net{i,1}(:,2);
InNN = 2*(SamIn'-repmat(MinSam,[1 num]))./repmat(MaxSam-MinSam,[1 num])-1;
%
G = Net{i,2}*InNN+repmat(Net{i,3},[1 num]);
U = tansig(G);
OutNN = Net{i,4}*U+Net{i,5};
%
SamOut = (OutNN+1)/2*(Net{i,6}(1,2)-Net{i,6}(1,1))+Net{i,6}(1,1);
SamOut = SamOut';
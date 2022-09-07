function output_vect=KeepOnly(vect,percentage)
% this function keeps only the highest 'percentage'% of the values in the
% input vecor
% 
% Inputs:  - vect: the input vector
%          - percentage: the desired percentage

%
% Outputs: -output_vect: output vector after keeping only the highest
%                        (percentage)% values, other values are set to zero
%
% This code was originally developped by Sahar Yassine
% contact: saharyassine94@gmail.com
%%
nToKeep=floor((percentage/100)*length(vect));
[sortedV,ind]=sort(vect,'descend');
indTokeep=ind(1:nToKeep);
output_vect=zeros(size(vect));
output_vect(indTokeep)=vect(indTokeep);

end
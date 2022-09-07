function [W,group]=SNF_Clust_With_hyperparams(features,K,alpha,nClust)
% This function performs the Similarity Network Fusion clustering using functions from the
% SNF matlab toolbox
%
% Inputs:  - features: input features of dimension [nSubj,nFeatures,nDataType]
%          - K: hyperparameter value of k-nearest neighbors 
%          - alpha: hyperparameter value of the scaling factor (mu) 
%          - nClust: the number of clusters
%
% Outputs: - W: fused similarity matrix of the pairwise hyperparameters
%               (k,alpha), dim [nSubj,nSub]
%          - group: final group affiliation of the subjects
%
% This code was originally developped by Sahar Yassine
% contact: saharyassine94@gmail.com

%%
T = 20; %Number of Iterations, usually (10~20)

[nSubj,~,nDataType]=size(features);

data_norm=zeros(size(features));
for i=1:nDataType
    data_norm(:,:,i)=Standard_Normalization(squeeze(features(:,:,i)));
end

%%%Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows;
%if the data is discrete, we recommend the users to use chi-square distance
Dist=zeros(nSubj,nSubj,nDataType);
for i=1:nDataType
    Dist(:,:,i)=dist2(squeeze(data_norm(:,:,i)),squeeze(data_norm(:,:,i)));
end



%%% construct similarity graphs
Ws_cell=cell(1,nDataType);
for i=1:nDataType
    Ws_cell{1,i}=affinityMatrix(squeeze(Dist(:,:,i)),K,alpha);
end

%Fuse the affinity matrices across dataTypes and apply spectral clustering
W = SNF(Ws_cell, K, T);
group = SpectralClustering(W,nClust);%%%the final subtypes information



end

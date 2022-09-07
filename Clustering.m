function [assignments,Ws_stable,group_affiliation,AssociationMat]=Clustering(features,K,alpha,clusters,display)
% This function performs the clustering pipeline using Similarity Network
% Fusion (SNF)
%
% Inputs:  - features: input features of dimension [nSubj,nFeatures,nDataType]
%          - K: hyperparameter vector for k-nearest neighbors 
%          - alpha: hyperparameter vector of the scaling factor (mu) 
%          - clusters: vector of number of cluster's solutions
%          - display: if 1 (default) ==> diplay is on, otherwise ==> no diplay
%
% Outputs: - assignments: group affiliation of the stable hyperparameters
%                         ,dim [nSubj, nStableOpt]
%          - Ws_stable : fused similarity matrices of the stable
%                        hyperparameters, dim [nSubj,nSub,nStableOpt]
%          - group_affiliation: final group affiliation of the subjects
%          - AssociationMat: co-assignment matrix fo the subjects
%                            indicating the probability of two patients being assigned to the
%                            same cluster across the stable clustering
%                            solutions, dim [nSubj,nSubj]
%          
%
% This code was originally developped by Sahar Yassine
% contact: saharyassine94@gmail.com

%%

if (nargin<5)
    display=1;
end

[nSubj,~,~]=size(features);
nOpt=length(K)*length(alpha);
nClusters=length(clusters);

Ks=zeros(nOpt,nClusters);
alphas=zeros(nOpt,nClusters);
groups=zeros(nSubj,nOpt,nClusters);
Ws=zeros(nSubj,nSubj,nOpt,nClusters);


for k=1:nClusters
    clust=clusters(k);
    count=1;
    for j=1:length(alpha)
        aa=alpha(j);
        for i=1:length(K)
            kk=K(i);
            Ks(count,k)=kk;
            alphas(count,k)=aa;
            [W,group]=SNF_Clust_With_hyperparams(features,kk,aa,clust);
            groups(:,count,k)=group;
            Ws(:,:,count,k)=W;
            count=count+1;
        end
    end
end

zrandsAvg=zeros(nOpt,nClusters);
zrandsThresh=zeros(nOpt,nClusters);
for i=1:nClusters
    in_ks=Ks(:,i);
    in_alphas=alphas(:,i);
    in_groups=groups(:,:,i);
    zrandsAvg(:,i)=Average_zrand(in_ks,in_alphas,in_groups);
    threshZ=KeepOnly(zrandsAvg(:,i),5);
    zrandsThresh(:,i)=threshZ;
end


assignments=[];
Ws_stable=[];
for j=1:nClusters
    for i=1:nOpt
        if(zrandsThresh(i,j)>0)
            group=groups(:,i,j);
            W=Ws(:,:,i,j);
            assignments=cat(2,assignments,group);
            Ws_stable=cat(3,Ws_stable,W);
        end
    end
end


[newM ,q, AssociationMat] = consensus_iterative(assignments');
group_affiliation=newM(1,:)';

if (display>0)
    displayClusters(AssociationMat,group_affiliation);
end


end
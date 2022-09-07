function accuracy=Robustness_of_Clusters(features,alpha,nclusters,clust_affiliation_orig)
% This function examines the robustness of the clusters after randomly
% removing 10% of the subjects, and comparing the obtained affiliation
% with original one using the Rand index
%
% Inputs: - features: input features of dimension [nSubj,nFeatures,nDataType]
%         - alpha: hyperparameter vector of the scaling factor (mu) 
%         - nclusters: number of desired clusters (number of clusters in
%                      the final group affiliation obtained from the
%                      clustering pipeline
%         - clust_affiliation_orig: the final group affiliation obtained from the
%                                   clustering pipeline
%
% Outputs: - accuracy : the final robustness percentage
%
% This code was originally developped by Sahar Yassine
% contact: saharyassine94@gmail.com
%%

nSubj=size(features,1);
nRemove=floor(0.1*nSubj);
accuracy=[];


for i=1:100
    mask=ones(nSubj,1);
    idRemove=randi([1,nSubj],nRemove,1);
    mask(idRemove)=0;
    [assignments,Ws_stable,fnewTemp,AssociationMat]=Clustering(features(mask>0,:,:),1:nSubj-nRemove-1,alpha,nclusters,0);
    ri=rand_index(clust_affiliation_orig(mask>0),fnewTemp);
    accuracy=[accuracy; ri];
    
end
accuracy =100 * mean(accuracy);

end
function zrandsAvg=Average_zrand(ks,alphas,groups)
%This function calculate the pairwise similarity index z-rand for each pair of (k,alpha) by
% averaging the z-rands between the clustering solution of the hyperparameter combination and 
% their four closest neighbours (k-1,alpha),(k+1,alpha),(k,alpha-stepAlpha),
% (k,alpha+stepAlpha)
%
% Inputs: - ks: vector of k used for all the optimization process,dim[nOpt,1]
%         - alphas: corresponding vector of alpha used for all the
%                   optimization process, dim[nOpt,1]
%         - groups: clustering solution for the corresponding pair of
%                   hyperparameters in (ks,alphas), dim [nubjects,nOpt]
%
% Outputs: -zrandsAvg: average z-rand for each pairwise od (k,alpha),
%                     dimension [nOpt,1]
%
% This code was originally developped by Sahar Yassine
% contact: saharyassine94@gmail.com

%%
nOpt=length(ks);
uniqueK=unique(ks);uniqueAlpha=unique(alphas);
sortedK=sort(uniqueK);sortedAlpha=sort(uniqueAlpha);
stepK=sortedK(2)-sortedK(1);stepAlpha=sortedAlpha(2)-sortedAlpha(1);
zrandsAvg=zeros(nOpt,1);

for i=1:nOpt
    k=ks(i);
    a=alphas(i);
    group=groups(:,i);
    nNeighbors=4;
    
    %1st neighbor
    k_n1=k; a_n1=a-stepAlpha;
    id_n1=find(ks==k_n1 & alphas==a_n1);
    if isempty(id_n1)
        zrand_n1=0;
        nNeighbors=nNeighbors-1;
    else
        group_n1=groups(:,id_n1);
        zrand_n1=zrand(group,group_n1);
    end
    
    
    %2nd neighbor
    k_n2=k; a_n2=a+stepAlpha;
    id_n2=find(ks==k_n2 & alphas==a_n2);
    if isempty(id_n2)
        zrand_n2=0;
        nNeighbors=nNeighbors-1;
    else
        group_n2=groups(:,id_n2);
        zrand_n2=zrand(group,group_n2);
    end
    
    
    %3rd neighbor
    k_n3=k-stepK; a_n3=a;
    id_n3=find(ks==k_n3 & alphas==a_n3);
    if isempty(id_n3)
        zrand_n3=0;
        nNeighbors=nNeighbors-1;
    else
        group_n3=groups(:,id_n3);
        zrand_n3=zrand(group,group_n3);
    end
    
    
    %4th neighbor
    k_n4=k+stepK; a_n4=a;
    id_n4=find(ks==k_n4 & alphas==a_n4);
    if isempty(id_n4)
        zrand_n4=0;
        nNeighbors=nNeighbors-1;
    else
        group_n4=groups(:,id_n4);
        zrand_n4=zrand(group,group_n4);
    end
    
    
    % average on neighbors
    zrandsAvg(i)=(zrand_n1+zrand_n2+zrand_n3+zrand_n4)/nNeighbors;
    
end



end
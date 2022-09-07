function [featuresPN,pvalsAll,significance_Matrix]=Discriminative_RSN_features(features,group_affiliation,RSN,cov,continuous_variables,pThresh)
% this function calculate the average relative power per resting
% state network and compute n-way anovas (covariate control) to retreive
% the significant network per dataType(frequency band) and their
% corresponding p-values, paiwise of networks/band who survive the
% Bonferroni correction are retained.
% 
% Inputs:  - features: input features of dimension [nSubj,nFeatures,nDataType]
%          - group_affiliation: the final group affiliation obtained from the
%                               clustering pipeline
%          - RSN: table containing the four resting state networks and the
%                 indexes of the corresponding brain regions (this table is
%                 saved as a .mat variable in the repositary folder
%          - cov: matrix containing the covariable to include in the anova
%                 analysis, dim [nSubj,nCov]
%          - continuous_variables: vector containing the indexes of the
%                                  continuous covariables
%          - pThresh: threshold p-value for significance (before bonferroni
%                     correction)
%
% Outputs: -featuresPN: structure of the mean features per resting state network for each subject
%          -pvalsAll: p-values resulted from the Anova analysis, dim[nNetwork, nDataType]
%          -significance_Matrix: fignificance matrix resulted from the
%                                Anova analysis, , dim[nNetwork, nDataType]
%
%          
%
% This code was originally developped by Sahar Yassine
% contact: saharyassine94@gmail.com
%%

nNetworks=size(RSN,2);
nDataType=size(features,3); %here represent the number of frequency bands
pvalsAll=zeros(nNetworks,nDataType);
significance_Matrix=zeros(nNetworks,nDataType);
nComparisons=nNetworks*nDataType;
pThresh_Bonf=pThresh/nComparisons;

% here covariables are as {age,Sex,Edu,LEDD}
group_cov=[group_affiliation,cov]; % the continuous variables are thus 2,4 and 5
continuous_variables=continuous_variables+1; % to change if needed

%Calculate the mean power per network
featuresPN.SMN=squeeze(mean(features(:,RSN.SMN>0,:),2));
featuresPN.DMN=squeeze(mean(features(:,RSN.DMN>0,:),2));
featuresPN.FTN=squeeze(mean(features(:,RSN.FTN>0,:),2));
featuresPN.FPN=squeeze(mean(features(:,RSN.FPN>0,:),2));



% n-way anovas (cov control) on SMN
for i=1:nDataType
     pvs=anovan(featuresPN.SMN(:,i),group_cov,'continuous',continuous_variables,'display','off');
     pvalsAll(1,i)=pvs(1);
     %check for significance
     if (pvs(1)<pThresh_Bonf)
         if(sum(pvs(2:end)<pThresh)==0)
             significance_Matrix(1,i)=1;
         end
     end
end

% n-way anovas (cov control) on DMN
for i=1:nDataType
     pvs=anovan(featuresPN.DMN(:,i),group_cov,'continuous',continuous_variables,'display','off');
     pvalsAll(2,i)=pvs(1);
     %check for significance
     if (pvs(1)<pThresh_Bonf)
         if(sum(pvs(2:end)<pThresh)==0)
             significance_Matrix(2,i)=1;
         end
     end
end

% n-way anovas (cov control) on FTN
for i=1:nDataType
     pvs=anovan(featuresPN.FTN(:,i),group_cov,'continuous',continuous_variables,'display','off');
     pvalsAll(3,i)=pvs(1);
     %check for significance
     if (pvs(1)<pThresh_Bonf)
         if(sum(pvs(2:end)<pThresh)==0)
             significance_Matrix(3,i)=1;
         end
     end
end

% n-way anovas (cov control) on FPN
for i=1:nDataType
     pvs=anovan(featuresPN.FPN(:,i),group_cov,'continuous',continuous_variables,'display','off');
     pvalsAll(4,i)=pvs(1);
     %check for significance
     if (pvs(1)<pThresh_Bonf)
         if(sum(pvs(2:end)<pThresh)==0)
             significance_Matrix(4,i)=1;
         end
     end
end


end

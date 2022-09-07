function finalEmbedding=Compute_embeddings(Ws_stable)
% This function performs the dimensionality reduction pipeline over the
% stable affinity matrices obtained from the clustering pipeline, using
% functions from the Brainspace embedding toolbox
%
% Inputs:  - Ws_stable : fused similarity matrices of the stable
%                        hyperparameters, dim [nSubj,nSub,nStableOpt]
%
% Outputs: - finalEmbedding : the final set of low-dimensional
%                             representation of the subjects across 5 final components
%          
%
% This code was originally developped by Sahar Yassine
% contact: saharyassine94@gmail.com

size_stable_W=size(Ws_stable,3);
grad=cell(size_stable_W,1);

for i=1:size_stable_W
    w=squeeze(Ws_stable(:,:,i));
    [embedding] = diffusion_mapping(w,5,0,0);
    grad{i,1}=embedding;
end

% align using Procrustes analysis
[aligned] = procrustes_alignment(grad);

%average aligned embedding 
finalEmbedding=[];
for i=1:length(aligned)
finalEmbedding=cat(3,finalEmbedding,aligned{1,i});
end
finalEmbedding=mean(finalEmbedding,3);

end
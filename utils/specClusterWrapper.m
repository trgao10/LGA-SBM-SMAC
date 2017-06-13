function [clusterLabel,embedPts] = specClusterWrapper(rIdx, cIdx, weightedAdj, numClusters, numKmeans, bandwidth, type)
%SPECTRALCLUSTERINGWRAPPER run spectral clustering multiple times and pick
%                          the best result
%   
% INPUTS:
%   weightedAdj ---------------------------- weighted adjacency matrix
%                                            (must be symmetric)
%   numClusters ---------------------------- number of clusters
%   numKmeans ------------------------------ number of runs
%   bandwidth ------------------------------ diffusion bandwidth paramerter
%                                            for the anisotropic kernel
%                                            (not useful if type == 'sim')
%   type ----------------------------------- indicating whether weightedAdj
%                                            contains similarity ('sim') or
%                                            dissimilarity ('dis') scores
%   debugFlag ------------------------------ print to screen
%
% OUTPUTS:
%   clusterLabel --------------------------- cell array of length
%                                            numClusters
%   
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 25, 2016
% 

vals = zeros(size(rIdx));
for j=1:length(vals)
    vals(j) = weightedAdj(rIdx(j), cIdx(j));
end

if strcmpi(type, 'dis')
    if isnumeric(bandwidth)
        simScore = exp(-vals/bandwidth);
    elseif ischar(bandwidth)
        simScore = exp(-vals/mean(vals));
    end
elseif strcmpi(type, 'sim')
    if (max(vals)-min(vals)) > 0
        simScore = (vals-min(vals))/(max(vals)-min(vals));
    else
        simScore = vals / max(vals);
    end
else
    error('unknown type');
end

W = sparse(rIdx,cIdx,simScore,size(weightedAdj,1),size(weightedAdj,2));
W = W+W';

Dvec = sum(W);
L = diag(1./sqrt(Dvec))*W*diag(1./sqrt(Dvec));
L = (L+L')/2;

[evecs, evals] = eig(L);
devals = diag(evals);
devals(devals < 0) = 0;
devals = 1 - devals;
[devals,evalidx] = sort(devals,'ascend');
evecs = evecs(:,evalidx);
embedPts = diag(1./sqrt(Dvec))*evecs(:,1:(numClusters+1))*sqrt(diag(1-devals(1:(numClusters+1))));
% embedPts = diag(1./sqrt(Dvec))*evecs(:,2:(numClusters+1))*sqrt(diag(1-devals(2:(numClusters+1))));

% embedPts = sign(diag(1./sqrt(Dvec))*evecs(:,2)*diag(sqrt(devals(2))));
% embedPts = diag(1./sqrt(Dvec))*evecs(:,2:(numClusters+1));
% embedPtsFull = diag(1./sqrt(Dvec))*evecs;
% embedPtsFull = diag(1./sqrt(Dvec))*evecs*sqrt(diag(1-devals));
% embedPts = diag(1./sqrt(Dvec))*evecs(:,2:(numClusters+1))*diag(sqrt(devals(2:(numClusters+1))));
% embedPts = diag(1./sqrt(Dvec))*evecs(:,1:(numClusters+1))*diag(sqrt(devals(1:(numClusters+1))));

cluster_idx = kmeans(embedPts,numClusters,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
clusterLabel = cell(1,numClusters);
for j=1:numClusters
    clusterLabel{j} = find(cluster_idx == j);
end

end

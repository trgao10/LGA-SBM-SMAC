%% setup parameters
clearvars
close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

syncableFlag = false;
debugFlag    = true;
newCaseFlag  = true;

%%% load a saved random instance for reproducibilty
if ~newCaseFlag
    casePath = 'reproduce/successCase04.mat';
end

NVec = [25,25];  %%% number of vertices in each cluster of the SBM ---- length of
                %%% this vector indicates the number of clusters
d = 3;   %%% dimension of the orthogonal group

%%% probabilities (p,q) for generating SBM graph
p = 0.5; %0.6;   %%% in-cluster connection probability
q = 0.4;   %% between-cluster connection probability

%%% probabilities (pp,qq) for adding noise to each block
pp = 1; %%% in-cluster noise probability will be 1-pp
qq = 0; %%% between-cluster noise probability will be 1-qq

%%% parameters controlling the behavior of LGA/SMAC
ROffType = 'SO';
maxIter = 10;
tol = 1e-8;
numKmeans = 200;
bandwidth = 1;
adjType = 'dis';
colorList = {'r','b','k','m'};
hsv = rgb2hsv(parula);
close(gcf);

%%%%%%% Spectral Relaxation
syncRoutine = @syncSpecRelax;
%%%%%%% SDP Relaxation: less scalable option
%%%%%%%%%%% if using SDP, must make sure CVX is in the path %%%%%%%%%%%%%%%
% syncRoutine = @syncSDPRelax;
% run ~/Documents/MATLAB/cvx/cvx_startup.m

numClusters = length(NVec);
params = struct('debugFlag', debugFlag,...
    'd', d,...
    'numClusters', numClusters,...
    'tol', tol,...
    'maxIter', maxIter,...
    'numKmeans', numKmeans,...
    'bandwidth', bandwidth,...
    'adjType', adjType,...
    'hsv', hsv,...
    'syncRoutine', syncRoutine,...
    'ROffType', ROffType);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% generate random graph and edge potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if newCaseFlag
    G = genSBMExpander(NVec, p, q, debugFlag);
    disp('Press any key to continue......');
    pause();
    
    % generate "ground truth" vertex and edge potentials, and contaminate the
    % edge potential on cross-cluster links with random orthonormal matrices
    vertPotCell = cell(1,sum(NVec));
    for j=1:sum(NVec)
        R = genRandGroupElement(ROffType,d);
        vertPotCell{j} = R;
    end
    
    edgePotCell = cell(size(G.adjMat));
    [rIdx, cIdx] = find(G.adjMat);
    for j=1:length(rIdx)
        edgePotCell{rIdx(j),cIdx(j)} = vertPotCell{rIdx(j)}*vertPotCell{cIdx(j)}';
    end
    
    %%% add noise to the edge potentials
    if ~syncableFlag
        for j=1:length(G.ccRowIdx)
            if rand(1) < 1-qq %%% perturb cross-cluster edge potentials with probability 1-qq
                R = genRandGroupElement(ROffType,d);
                edgePotCell{G.ccRowIdx(j),G.ccColIdx(j)} = vertPotCell{G.ccRowIdx(j)}*R'*vertPotCell{G.ccColIdx(j)}';
                edgePotCell{G.ccColIdx(j),G.ccRowIdx(j)} = edgePotCell{G.ccRowIdx(j),G.ccColIdx(j)}';
            end
        end
        
        for j=1:length(G.icRowIdx)
            if rand(1) < 1-pp %%% perturb in-cluster edge potentials with probability 1-pp
                R = genRandGroupElement(ROffType,d);
                edgePotCell{G.icRowIdx(j),G.icColIdx(j)} = vertPotCell{G.icRowIdx(j)}*R'*vertPotCell{G.icColIdx(j)}';
                edgePotCell{G.icColIdx(j),G.icRowIdx(j)} = edgePotCell{G.icRowIdx(j),G.icColIdx(j)}';
            end
        end
    else
        R = genRandGroupElement(ROffType,d);
        for j=1:length(G.ccRowIdx)
            edgePotCell{G.ccRowIdx(j),G.ccColIdx(j)} = vertPotCell{G.ccRowIdx(j)}*R'*vertPotCell{G.ccColIdx(j)}';
            edgePotCell{G.ccColIdx(j),G.ccRowIdx(j)} = edgePotCell{G.ccRowIdx(j),G.ccColIdx(j)}';
        end
    end
else
    %%%%% load previous case
    if exist('casePath', 'var')
        load(casePath);
        params.debugFlag = debugFlag;
    else
        error('casePath unspecified');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% run SynCut and collect result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debugFlag
    [params.colorList] = deal(colorList);
    [params.vertPotCell] = deal(vertPotCell);
end

%%%% Can use SynCut instead of SynCut, but SynCut does not provide
%%%% visualization functionality for debugging.
rslt = SMAC_TwoWay(G, edgePotCell, params);

%%%% very ad-hoc error counts and error rates
%%%% only works for two clusters!
rslt.errCounts = min(sum(rslt.clusterLabel{1} > NVec(1))+sum(rslt.clusterLabel{2} <= NVec(1)),...
    sum(rslt.clusterLabel{1} <= NVec(1))+sum(rslt.clusterLabel{2} > NVec(1)));
rslt.errRate = rslt.errCounts / size(G.adjMat, 1);

%%%% save this run for debugging
save('reproduce/currCase.mat', 'G', 'vertPotCell', 'edgePotCell', 'ROffType');



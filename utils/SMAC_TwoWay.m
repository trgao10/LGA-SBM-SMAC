function [rslt] = SMAC_TwoWay(G, edgePotCell, params)
%SMAC_TWOWAY
%   INPUTS:
%     G ---------------------------- graph
%     edgePotCell ------------------ edge potential in cell array format
%     params ----------------------- struct holding all parameters
%
%   OUTPUTS:
%     rslt ------------------------- struct holding all outputs
%
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: June 12, 2017
%

debugFlag = getoptions(params, 'debugFlag', false);
d = getoptions(params, 'd', Inf);
numClusters = getoptions(params, 'numClusters', 2);
numKmeans = getoptions(params, 'numKmeans', 200);
bandwidth = getoptions(params, 'bandwidth', 1);
adjType = getoptions(params, 'adjType', 'dis');
syncRoutine = getoptions(params, 'syncRoutine', @syncSpecRelax);
ROffType = getoptions(params, 'ROffType', 'O');

if numClusters > 2
    warning('For numClusters > 2, use the routine "SynCut" instead.');
    error('This routine only demos the details of SynCut using K=2.');
end

if debugFlag    
    try
        vertPotCell = params.vertPotCell;
    catch
        error('with debugFlag=true, must provide ground truth vertPotCell');
    end
    try
        colorList = getoptions(params, 'colorList', {'r','b','k','m'});
    catch
        error('with debugFlag=true, must provide ground truth colorList');
    end
    
    [unGCL, GCL, GCL_Dvec, GCL_W] = assembleGCL(G.adjMat, edgePotCell, d);
    
    figure('Position',[30,550,560,420]);
    plot(graph(G.adjMat), 'XData', G.V(:,1), 'YData', G.V(:,2));
    axis equal
    axis([0,numClusters,0,1]);
    hold on
    for j=1:length(G.ccRowIdx)
        staPtCoords = G.V(G.ccRowIdx(j),:);
        endPtCoords = G.V(G.ccColIdx(j),:);
        line([staPtCoords(:,1);endPtCoords(:,1)],...
            [staPtCoords(:,2);endPtCoords(:,2)],'Color','r');
    end
    title(sprintf('Spectral Gap = %.2f, CCE/TTE = %d/%d',...
        G.specGap, length(G.ccRowIdx), sum(G.adjMat(:))/2),'Interpreter','latex');
    
    vertPotMat = cat(1, vertPotCell{:});
    [GroundTruthPerEdgeFrustVec, GroundTruthPerEdgeFrustMat] =...
        getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, vertPotCell);
    rescaled_vertPotMat = diag(sqrt(GCL_Dvec))*vertPotMat;
    fprintf('[GroundTruth] Rayleigh quotient (normalized GCL) = %f\n',...
        trace(rescaled_vertPotMat'*GCL*rescaled_vertPotMat));
    fprintf('[GroundTruth] Rayleigh quotient (unnormalized GCL) = %f\n',...
        trace(vertPotMat'*unGCL*vertPotMat));
    fprintf('[GroundTruth] total frustration = %f\n',...
        sum(GroundTruthPerEdgeFrustVec));
end

wAdjMat = G.adjMat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 1. Synchronization (global)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unwGCL, wGCL, wGCL_Dvec] = assembleGCL(wAdjMat, edgePotCell, d);
[RelaxSolCell, RelaxSolMat] = syncRoutine(unwGCL, d, wGCL_Dvec, ROffType);
[RelaxSolPerEdgeFrustVec, RelaxSolPerEdgeFrustMat] =...
    getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, RelaxSolCell);

if debugFlag
    rescaled_RelaxSolMat = diag(sqrt(GCL_Dvec))*RelaxSolMat;
    fprintf('[RelaxSol] Rayleigh quotient (normalized GCL) = %f\n',...
        trace(rescaled_RelaxSolMat'*wGCL*rescaled_RelaxSolMat));
    fprintf('[RelaxSol] Rayleigh quotient (unnormalized GCL) = %f\n',...
        trace(RelaxSolMat'*unwGCL*RelaxSolMat));
    fprintf('[RelaxSol] total frustration = %f\n',...
        sum(RelaxSolPerEdgeFrustVec));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 2. Spectral Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rIdx, cIdx] = find(G.adjMat);
while true
    [clusterLabel, embedPts] =...
        specClusterWrapper(rIdx, cIdx, RelaxSolPerEdgeFrustMat,numClusters,...
        numKmeans,bandwidth,adjType);
    %%%%% we expect spectral clustering to partition the graph into two
    %%%%% connected components
    connTest = zeros(size(clusterLabel));
    numConnComp = zeros(size(clusterLabel));
    for j=1:length(connTest)
        connTest(j) = ~isempty(find(sum(G.adjMat(clusterLabel{j},clusterLabel{j})) == 0, 1));
        numConnComp(j) = length(unique(conncomp(graph(G.adjMat(clusterLabel{j},clusterLabel{j})))));
    end
    if (sum(connTest) == 0) && (all(numConnComp == 1))
        break
    end
    error('error: spectral clustering does not produce connected components');
end

if debugFlag
    %%%% for more consistent visualization effect, always set left
    %%%% cluster to be red an right cluster to be blue
    %%%% Note this will only work for K=2 (binary clustering)!
    NVec = G.NVec;
    if sum(clusterLabel{1} <= NVec(1)) < sum(clusterLabel{2} <= NVec(1))
        tmp = clusterLabel{1};
        clusterLabel{1} = clusterLabel{2};
        clusterLabel{2} = tmp;
        clear tmp
    end
    
    %%%% plot spectral clustering results
    if ~exist('specClusteringFigure', 'var')
        specClusteringFigure = figure('Position',[1200,550,1000,400]);
    else
        figure(specClusteringFigure);
    end
    
    for j=1:length(clusterLabel)
        scatter(G.V(clusterLabel{j},1),G.V(clusterLabel{j},2),...
            20,colorList{j},'filled');
        if j==1
            hold on
        end
    end
    axis equal
    axis([0,numClusters,0,1]);
    line([1,1],[0,1],'Color','g');
end
    
rslt = struct('G', G, 'params', params);
[rslt.edgePotCell] = deal(edgePotCell);
[rslt.clusterLabel] = deal(clusterLabel);
rslt.embedPts = embedPts;

end

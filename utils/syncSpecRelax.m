function [solCell, solMat] = syncSpecRelax(GCL, d, GCL_Dvec, ROffType)
%SYNCSPECRELAX Solve synchronization problem using spectral relaxation
% Based on the work [BSS2013]:
% =========================================================================
%   Bandeira, Afonso S., Amit Singer, and Daniel A. Spielman.
%   "A Cheeger Inequality for the Graph Connection Laplacian."
%   SIAM Journal on Matrix Analysis and Applications 34.4 (2013):1611-1630.
% =========================================================================
%
% INPUTS:
%   GCL ----------------------- Graph Connection Laplacian
%   d   ----------------------- (fixed) dimension of each block
%   GCL_Dvec ------------------ vector encoding the D matrix in the formula
%                               GCL=I-D^{-1/2}*W*D^{-1/2}
%   ROffType ------------------ type of round-off to be performed
%                               [ 'O' ] | 'SO' | 'Perm'
%
% OUTPUTS:
%   solCell ------------------- solution in cell format (cell array of
%                               length n)
%   solMat -------------------- solution in matrix format (matrix of size
%                               (n x numClusters)-by-n
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: June 11, 2017
%

if nargin < 4
    ROffType = 'O';
    disp('ROffType not specified.... set to "O" by default');
end

disp(['ROffType = ' ROffType]);

[evecs,evals] = eig((GCL+GCL')/2);
evals = diag(evals);
[~,idx] = sort(evals);
evecs = evecs(:,idx);
relaxSol = evecs(:,1:d);
% rescaled_relaxSol = diag(1./sqrt(GCL_Dvec))*relaxSol;
rescaled_relaxSol = sqrt(size(GCL,1)/d) * relaxSol;

fprintf('++++++++++++++++++++\n');
fprintf('relaxed frustration before rounding off: %f\n', trace(rescaled_relaxSol'*GCL*rescaled_relaxSol));

solCell = cell(1,length(GCL_Dvec)/d);
solMat = zeros(size(rescaled_relaxSol));
for j=1:length(GCL_Dvec)/d
    tmpIdx = ((j-1)*d+1):(j*d);
    tmpBlock = rescaled_relaxSol(tmpIdx,:);
    if strcmpi(ROffType, 'O')
        [U,~,V] = svd(tmpBlock);
        solCell{j} = U*V';
    elseif strcmpi(ROffType, 'SO')
        [U,~,V] = svd(tmpBlock);
        solCell{j} = U*[[eye(d-1),zeros(d-1,1)];zeros(1,d-1),sign(det(U*V'))]*V';
    elseif strcmpi(ROffType, 'Perm')
        tmpP = assignmentoptimal(tmpBlock-min(tmpBlock(:)));
        solCell{j} = full(sparse(1:length(tmpP), tmpP, ones(1,length(tmpP))));
    else
        error('unknown ROffType');
    end
    solMat(tmpIdx,:) = solCell{j};
end

fprintf('relaxed frustration after rounding off: %f\n', trace(solMat'*GCL*solMat));
fprintf('++++++++++++++++++++\n');

end


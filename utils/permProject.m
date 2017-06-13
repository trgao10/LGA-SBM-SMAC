function [ pPerm, objVal ] = permProject( x, Aeq, beq, partialD )
%===============================================================
% module:
% ------
% permProject.m
%
% project:
% -------
% Thesis: Algorithm for solving trust region fashion
%
% Description:
% -----------
% This function finds the closet permutation to a given matrix, i.e solves
% the LP:
%   max x^Tp
%   s.t A_{eq}p=b_{eq}
%   0<= p_i <= 1
%
% Inputs:
% ------
% x: input matrix as a vector column stuck
% Aeq: matrix as in description, describes 1^TX=1^T, X1=1.
% beq: vector as in description, describes 1^TX=1^T, X1=1.
%
% Outputs:
% -------
% pPerm: the closet permutation to a given matrix
% objVal: objective value
%
% Authors: Itay kezurer
% Copyright 2014 Kezurer
% Revision: 1.0 Date: 12/10/2014 21:29
%===============================================================
%--------------------------------------------
% Initialization
%--------------------------------------------
consTerm = 10^(-5);
options.Display = 'off';
% randVec = consTerm * rand(size(x));
%============================================
%--------------------------------------------
% LP solver
%--------------------------------------------
% with minus sign because in matlab it's min
% [pPerm, objVal] = linprog(-(x + randVec),[],[],Aeq,beq,zeros(size(x)),ones(size(x)),[],options);
if ~exist('partialD','var') || (partialD==sqrt(length(x)))
    % standard projection on permutations
    [pPerm, objVal] = bintprog(-x,[],[],Aeq,beq,[],options);
else
    % projection on partial permutations
    [pPerm, objVal] = bintprog(-x,Aeq,beq,ones(1,length(x)),partialD,[],options);
end
% [pPerm, objVal] = linprog(-x,[],[],Aeq,beq,zeros(size(x)),ones(size(x)),[],options);
% pPerm_prev = pPerm;
% pPerm = round(pPerm);
% if norm(pPerm_prev(:)-pPerm(:),'inf')>consTerm
%     warning('### some rounding problem has occured in permProject ###');
% end
%============================================

end


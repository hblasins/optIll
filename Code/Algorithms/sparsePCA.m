function [X, wghts, S, hist] = sparsePCA(S, d, lambda, epsilon, basis, varargin)

% [X, wghts, S, hist] = sparsePCA(S, d, lambda, epsilon, bases, ...)
%
% This function implements the complete sparse nonnegative PCA algorithm
% (algorthm 2 in the manuscript). The user specified the number of
% projection directions to compute. At each step the optimal sparse,
% nonnegative direction is computed using the nonnegative Fantope PCA
% algorithm. Next the covariance matrix is deflated with the Mackey method.
%
% Inputs:
%    S - the covariance matrix of the set of reflectance or 
%    reflectance-responsivity spectra.
%    d - the number of projection directions to compute.
%    lambda - a scalar controling the strength of the sparsity penalty.
%    epsilon - a scalar controling the desired accuracy of the ADMM
%    solution.
%    basis - an n x k matrix of k spectral basis functions. The optimal illuminant
%    lies in the span of columns of bases. 
%
% Outputs:
%    X - a n x d matrix of sparse, nonnegative PCA directions.
%    wghts - a k x d matrix of weights such that X=bases*weights.
%    S - a n x n covariance matrix after d deflation steps (i.e. part of
%    the covariance matrix 'unaccounted' by the d nonnegative sparse PCA
%    projections).
%    hist - a d element cell array containing the ADMM convergence information 
%    (primal and dual residuals) for every nonegative sparse PCA direction. 
%
% Copyright, Henryk Blasinski 2017


n = size(S,1);
B = eye(n);
k = size(basis,2);

wghts = zeros(k,d);
X = zeros(n,d);
hist = cell(d,1);



for i=1:d

    % Find a rank-1 matrix, then use eigendecomposition to find the optimal
    % projection vector.
    [R, ~, hist{i}] = nonnegFantopePCA( S, 1, lambda, epsilon, basis, varargin{:} );
    [Pc, eigVals] = eig(R);
    [~, indx] = sort(diag(eigVals),'descend');
    Vhat = Pc(:,indx(1));
    
    % We assume that very small values are effectivel zero.
    Vhat(abs(Vhat) < 1e-4) = 0;

    qt = B*Vhat;
    Pihat = qt*qt';
    S = (eye(n) - Pihat)*S*(eye(n) - Pihat);
    B = B*(eye(n) - Pihat);
    
    tmp = Vhat;
    
    % Finally since xx^T = (-x)(-x^T) we need to make sure that we have a
    % vector with all nonnegative rather than nonpositive entries.
    sign = 1;
    if max(tmp) < max(-tmp) 
        sign = -1;
    end
    
    X(:,i) = sign*Vhat;
    
    w = basis\(sign*Vhat);
    wghts(:,i) = w;
    
end

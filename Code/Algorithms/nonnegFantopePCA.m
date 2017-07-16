function [ X, W, hist ] = nonnegFantopePCA( S, d, lambda, epsilon, B, varargin )

% [ X, W, hist ] = nonnegFantopePCA( S, d, lambda, epsilon, B, ... )
%
% This function computes the sparse nonnegative PCA directon using the ADMM
% algorithm described in the manuscript (Algorithm 1).
%
% Inputs:
%    S - the covariance matrix of the set of reflectance or 
%    reflectance-responsivity spectra.
%    d - the rank of the solution X. If d=1 then X is rank-1 and is therefore an outer
%    product xx^T, where x is the optimal, nonegative sparse PCA direction.
%    lambda - a scalar controling the strength of the sparsity penalty.
%    epsilon - a scalar controling the desired accuracy of the ADMM
%    solution.
%    B - an n x k matrix of k spectral basis functions. The optimal illuminant
%    lies in the span of columns of B. 
%
% Outputs:
%    X - a n x n, rank d (approximately) matrix
%    W - a k x k matrix representing X in terms of basis function B weights. 
%    hist - a 'history' structure that contains the values of the primal
%    and dual residuals, rho and the number of preconditioned conjugate gradient
%    descent steps required to solve the least-squares problem. All these
%    parameters are stored for each ADMM algorithm step.
%
% Copyright, Henryk Blasinski 2017.

p = inputParser;
p.KeepUnmatched = true;
p.addOptional('maxIter',10000);
p.addOptional('mu',10);
p.addOptional('tauIncr',2);
p.addOptional('tauDecr',2);
p.addOptional('adjustRho',true);
p.addOptional('Verbose',0);
p.parse(varargin{:});
inputs = p.Results;



% These variables are associated with Fantope projection penalty
Y1 = zeros(size(B,1));
Y1minus = zeros(size(B,1));
U1 = zeros(size(B,1));

% These with nonnegativity and sparsity
Y2 = zeros(size(B,1));
Y2minus = zeros(size(B,1));
U2 = zeros(size(B,1));


W = zeros(size(Y1));

%%
BB = kron(B,B);

hist.rho = ones(inputs.maxIter+1,1);
hist.prRes = zeros(inputs.maxIter,1);
hist.dualRes = zeros(inputs.maxIter,1);
hist.pcgIter = zeros(inputs.maxIter,1);

fprintf('Running non-negative sparse PCA arbitrary ...');
tic
for i=1:inputs.maxIter
    
    % Optimize over W
    % PCG approach (doesn't work if B is fat)
    %{
    b = applyAtb(B,Y1-U1,Y2-U2);
    AtAhndl = @(x) applyAtA(x,B);
    
    tic;
    [W, ~, relres, hist.pcgIter(i)] = pcg(AtAhndl,b,1e-10,10000,[],[],W(:));
    t1 = toc;
    
    W = reshape(W,sqrt(size(W,1))*ones(1,2));
    %}
    b = (Y1-U1 + Y2-U2)/2;
    W = BB\b(:);
    W = reshape(W,sqrt(size(W,1))*ones(1,2));
    
    
    % First optimize over Fantope projection
    Y1 = real(fantopeProjection(B*W*B' + U1 + S/hist.rho(i),d));
        
    % Hinge threshold Y2 - sparsity and nonnegativity
    tmp = B*W*B'+U2;
    Y2 = (tmp > 0).*max((abs(tmp)-lambda/hist.rho(i)),0);
    
    % Update dual variables
    U1 = U1 + B*W*B' - Y1;
    U2 = U2 + B*W*B' - Y2;
    
    
    hist.prRes(i) = norm(B*W*B'-Y1,'fro')^2 + norm(B*W*B'-Y2,'fro')^2;
    hist.dualRes(i) = (hist.rho(i)^2)*(norm(Y1 - Y1minus,'fro')^2 + norm(Y2 - Y2minus,'fro')^2);
    
    
    if (hist.prRes(i) > inputs.mu*hist.dualRes(i)) && (inputs.adjustRho == true);
        hist.rho(i+1) = hist.rho(i)*inputs.tauIncr;
    else
        if (hist.dualRes(i) > inputs.mu*hist.prRes(i))
            hist.rho(i+1) = hist.rho(i)/inputs.tauDecr;
        else
            hist.rho(i+1) = hist.rho(i);
        end
    end
    
    % Now we need to update U because this is not the actual dual variable
    % but Z = rho*U is.
    U1 = U1*hist.rho(i)/hist.rho(i+1);
    U2 = U2*hist.rho(i)/hist.rho(i+1);

      
    if max(hist.prRes(i),hist.dualRes(i)) <= d*epsilon,
        break; 
    end;
    
    if inputs.Verbose >= 2
        fprintf('Iteration %i, condition %f\n',i,max(hist.prRes(i),hist.dualRes(i)));
    end
    % fprintf('     -> PCG (%f), err %f, nIter %i\n',t1,relres,hist.pcgIter(i));
    Y1minus = Y1;
    Y2minus = Y2;
     
end
X = Y1;
fprintf(' Done! (%i/%i iterations in %f sec, residual %f)\n',i,inputs.maxIter,toc,max(hist.prRes(i),hist.dualRes(i)));
hist.prRes = hist.prRes(1:i);
hist.dualRes = hist.dualRes(1:i);
hist.rho = hist.rho(1:i);

% figure; plot([prRes dualRes]);

end

function res = applyAtA(x,B)

    n = sqrt(length(x));
    x = reshape(x,[n n]);
    res = 2*(B'*B)*x*(B'*B);
    res = res(:);

end

function res = applyAtb(B,b1,b2)

    res = B'*(b1+b2)*B;
    res = res(:);
    
end

function P = fantopeProjection(X,d)

[V, D] = eig(X);
gamma = findTheta(diag(D),d);
P = V*diag(gamma)*V';

end


function [gamma, theta] = findTheta(eigVals,d)

u = max(eigVals);
l = min(eigVals) - 1;

for i=1:1000
    theta = 0.5*(u+l);
    gamma = min(max(eigVals - theta,0),1);
    fVal = sum(gamma);    
    
    if fVal < d
        u = theta;
    else
        l = theta;
    end
    
    if  abs(fVal - d) <= 1e-6, break; end;
end
end
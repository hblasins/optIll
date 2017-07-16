function [ phi, predLabels, W, b, hist ] = multiclassSVM( data, labels, basis, varargin)

% [ phi, predLabels, W, b, hist ] = multiclassSVM( data, labels, basis, ...)
%
% This function implements the supervised optimal illuminant selection
% algorithm as described in section 4 (Supervised framework) of the
% manuscript.
% 
% Briefly, the approach incorporates linear image formation model and
% optimizes jointly over linear SVM classifier weights and optimal
% illuminants.
%
% Note: This code uses CVX to solve convex optimization problems, it will
% not scale well with the problem size.
%
% Inputs:
%
%    data - a n x p array of p different spectral quantities.
%    labels - a p element vector of labels assigned to each spectrum.
%    basis - an n x k matrix of k basis functions. The computed optimal
%    illuminants will lie in the span of the columns of basis.
%
% Outputs:
%    phi - a n x nChannels matrix contaning optimal illuminants.
%    predLabels - the labels assigned to each of the p input spectra at the
%    end of training.
%    W - a nChannels x nClasses matrix containing linear SVM classifier
%    boundary hyperplane parameters.
%    b - a nClasses x 1 vector contaning linear SVM classifier boundary offsets.
%    hist - a structure containing information about convergence of the
%    biconvex optimization algorithm.
%
% Copyright, Henryk Blasinski 2017

p = inputParser;
p.KeepUnmatched = true;
p.addOptional('multiclassMaxIter',20);
p.addOptional('nChannels',3);
p.addOptional('multiclassC',1);
p.addOptional('alpha',0.0);
p.addOptional('basisWeights',[]);
p.addOptional('eps',1e-6);
p.addOptional('multiclassPenalty','sparsity');

p.parse(varargin{:});
inputs = p.Results;


classLabels = unique(labels);
k = length(classLabels);
l = size(data,2);
n = size(data,1);


% We have to initialize the basis Weights with a feasible (i.e. nonnegative
% responsivity curves.

basisWeights = randn([size(basis,2),inputs.nChannels]);


cvx_begin quiet
    variable t(size(basis,2),inputs.nChannels)
    minimize norm(t(:) - basisWeights(:))
    subject to
        basis*t >= 0
cvx_end

basisWeights = t;



basisWeights = basisWeights*diag(1./sum(basis*basisWeights));


hist.basisWeightsInit = basisWeights;
hist.optVal = [];
hist.acc = [];

for j=1:inputs.multiclassMaxIter
    
    phi = basis*basisWeights;
    b = zeros(k,1);
    
    % Optimize over SVM weights
    cvx_begin quiet
        variables W(inputs.nChannels,k) ksi(k-1,l) b(k,1)
        switch inputs.multiclassPenalty
            case 'sparsity'
                minimize 0.5*sum(norms(W,2,1)) + inputs.multiclassC*sum(sum(ksi)) + inputs.alpha*sum(norms(basis*basisWeights,1,1))
            case 'smoothness'
                R = [eye(n-1) zeros(n-1,1)] - [zeros(n-1,1) eye(n-1)];
                minimize 0.5*sum(norms(W,2,1)) + inputs.multiclassC*sum(sum(ksi)) + inputs.alpha*sum(norms(R*basis*basisWeights,2,1))
        end

        subject to
            ksi >= 0;
            for i=1:l
                Wbar = W;
                Wbar(:,labels(i)) = [];
                bbar = b;
                bbar(labels(i)) = [];
                W(:,labels(i))'*repmat(phi'*data(:,i),[1 k-1]) + b(labels(i)) >= (phi'*data(:,i))'*Wbar + 2 - ksi(:,i)' + bbar';
            end
    cvx_end
    
    hist.optVal = [hist.optVal cvx_optval];
    
    [~, predSVM] = max((phi'*data)'*W + repmat(b',[l 1]),[],2);
    hist.acc = [hist.acc sum(predSVM' == labels)/length(labels)];
    fprintf('>> Iteration %i/%i:\n',j,inputs.multiclassMaxIter);
    fprintf('>>    SVM obj %.3f, accuracy: %.2f (%i/%i)\n',cvx_optval,hist.acc(end),sum(predSVM'==labels),length(labels));
    
    % Optimize over projection matrix phi
    cvx_begin quiet
        variables basisWeights(size(basis,2),inputs.nChannels) ksi(k-1,l) b(k,1)
        switch inputs.multiclassPenalty
            case 'sparsity'
                minimize 0.5*sum(norms(W,2,1)) + inputs.multiclassC*sum(sum(ksi)) + inputs.alpha*sum(norms(basis*basisWeights,1,1))
            case 'smoothness'
                R = [eye(n-1) zeros(n-1,1)] - [zeros(n-1,1) eye(n-1)];
                minimize 0.5*sum(norms(W,2,1)) + inputs.multiclassC*sum(sum(ksi)) + inputs.alpha*sum(norms(R*basis*basisWeights,2,1))
        end
        subject to
            ksi >= 0;
            for i=1:l
                Wbar = W;
                Wbar(:,labels(i)) = [];
                bbar = b;
                bbar(labels(i)) = [];
                W(:,labels(i))'*repmat((basis*basisWeights)'*data(:,i),[1 k-1]) + b(labels(i)) >= ((basis*basisWeights)'*data(:,i))'*Wbar + 2 - ksi(:,i)' + bbar';
            end
            basis*basisWeights >= 0;
            sum(basis*basisWeights) == 1;
    cvx_end
    
    hist.optVal = [hist.optVal cvx_optval];
    
    
    phi = basis*basisWeights;
    [~, predLabels] = max((phi'*data)'*W + repmat(b',[l 1]),[],2);
    
    hist.acc = [hist.acc sum(predLabels' == labels)/length(labels)];
    fprintf('>>    Cam obj %.3f, accuracy: %.2f (%i/%i)\n',cvx_optval, hist.acc(end),sum(predLabels'==labels),length(labels));
    
    if abs(hist.optVal(end-1) - hist.optVal(end)) < inputs.eps,
        break;
    end
    
end

predLabels = predLabels';

end


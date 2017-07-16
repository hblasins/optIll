% This a test script for a sparse nonnegative PCA problem as defined by
% (12) in the manuscript. We use the proposed nonnegative Fantope PCA projection and
% compare the result of the computation to that obtained With a general
% purpose solver such as cvx.
% 
% For tractability we consider small problem sizes only.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

%% Problem data 
rng(1);

alpha = 0.1;

n = 13; % #wavelengths
k = 5;  % #basis functions

refl = rand(100,n);
S = cov(refl);

B = rand(k,n)';

%% CVX implementation

cvx_begin
    variables Wcvx(k,k)
    maximize trace(S*B*Wcvx*B') - alpha*sum(norms(B*Wcvx*B',1))
    subject to
        B*Wcvx*B' == semidefinite(n)
        (eye(n) - B*Wcvx*B') == semidefinite(n)
        trace(B*Wcvx*B') == 1
        B*Wcvx*B' >= 0
cvx_end

Xcvx = B*Wcvx*B';

%% Nonegative Fantope implementation

[ Xfp, Wfp, hist ] = nonnegFantopePCA( S, 1, alpha, eps, B);

%% Scatter plots

figure;
hold on; grid on; box on;
plot([hist.prRes hist.dualRes]);
title('ADMM convergence');
xlabel('Iteration');
legend('pr. Res','dual Res');


figure;
hold on; grid on; box on;
plot(Xcvx(:),Xfp(:),'.','markerSize',10);
title('Scatter plot');
xlabel('CVX');
ylabel('ADMM');




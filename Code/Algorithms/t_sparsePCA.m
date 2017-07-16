% This a test script for the full nonnegative sparse PCA algorithm
% 
% For tractability we consider small problem sizes only.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

%% Problem data 
rng(1);

lambda = 0.1;
epsilon = eps;
d = 3;

n = 13; % #wavelengths
k = 5;  % #basis functions

refl = rand(100,n);
S = cov(refl);

B = rand(k,n)';

%% Nonnegative, sparse PCA

[X, wghts, S, hist] = sparsePCA(S, d, lambda, epsilon, B);



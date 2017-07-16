% This a test script for the supervised illuminant selection algorithm
% 
% For tractability we consider small problem sizes only.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

%% Problem data 
rng(1);

n = 13; % #wavelengths
k = 5;  % #basis functions

data = rand(n, 100);
labels = randi(3,[1 100]);

B = rand(k,n)';

%% Nonnegative, sparse PCA

[ phi, predLabels, W, b, hist ] = multiclassSVM( data, labels, B);



% Set up Matlab path to work with the optimal illuminant selection code.
%
% Copyright, Henryk Blasinski 2017
close all;
clear all;
clc;

[codePath, parentPath] = olRootPath;

addpath(codePath);

addpath(fullfile(codePath,'Algorithms'));
addpath(fullfile(codePath,'Utilities'));


if exist(fullfile(codePath,'Devices'),'dir')
    addpath(fullfile(codePath,'Devices'));
    addpath(fullfile(codePath,'Devices','Flea3'));
    addpath(fullfile(codePath,'Devices','LEDCube'));
    addpath(fullfile(codePath,'Devices','PR715'));
end

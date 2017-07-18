% This script uses the PR715 to measure the spectral power distributions of
% the individual LEDs in the LEDCube.
%
% It requires that the hardware (LEDCube and PR715) be connected to the
% computer.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

[~, parentPath] = olRootPath();

%% Measure spectra

nChannels = 11;
pr715Com = 4;
LEDCubeCom = 'com6';


pr = pr715init(pr715Com);
data = zeros(173,nChannels);

for i=1:nChannels
    receipe=zeros(11,1);
    receipe(i)=1023;
    TL_Cube_lightRecipe(LEDCubeCom,[0, 0, 0, 0]',receipe,2); 
    
    [data(:,i), wave] = pr715spectrum(pr);
end

fclose(pr);

figure;
plot(wave,data);
xlabel('Wavelength, nm');

%% Save
fName = fullfile(parentPath,'Parameters','LEDCubeSpectra.mat');
ieSaveSpectralFile(wave,data,'LEDCube spectra',fName);


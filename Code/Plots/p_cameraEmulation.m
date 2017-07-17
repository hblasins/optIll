% This script demonstrates the similarity between real and emulated RGB
% cameras. It produces Fig. 2 in the paper.
%
% Copyright, Henryk Blasinski 2017

close all;
clear all;
clc;

ieInit;

[codePath, parentPath] = olRootPath();
destDir = fullfile(parentPath,'TestFigures');
if isempty(destDir) == 0 && exist(destDir,'dir') == 0
    mkdir(destDir);
end

%% Load and crop the emulated image

fName = fullfile(parentPath,'Images','GoPro5_6500K.mat');
data = load(fName);

figure;
hold on; grid on; box on;
plot(data.cameraIll);
plot(data.ledCameraApprox,'--');

x = 420;
y = 967;
delta = 5;

% Perform white balancing using the white patch in the Macbeth chart.
wp = data.Img(y-delta:y+delta,x-delta:x+delta,:);
wp = reshape(wp,[(2*delta+1)^2, 3]);
wp = mean(wp);

imgVec = reshape(data.Img,[size(data.Img,1)*size(data.Img,2),3]);
imgCorrVec = imgVec*diag(1./wp(:));

imgCorr = reshape(imgCorrVec,[size(data.Img,1), size(data.Img,2),3]);
imgCorr = imgCorr/max(imgCorr(:));
figure; imshow(imgCorr.^(1/2.2));

simCrop = imgCorr(678:1005,380:886,:);
simCrop = 1.4*simCrop/max(simCrop(:));
figure; imshow(simCrop.^(1/2.2));

% Extract Macbeth patch intensities.
cp = [391 997;873 1006;873 697;398 695];

sensor = sensorCreate('monochrome');
sensor = sensorSet(sensor,'size',[size(data.Img,1) size(data.Img,2)]);
emulatedMacbethSamples = zeros(24,3);

for i=1:3
   sensor = sensorSet(sensor,'volts',imgCorr(:,:,i));
   ieAddObject(sensor);
   sensorWindow();
   
   [samples, ~, ~, cp] = chartSelect(sensor,1,1,4,6,cp); 
   emulatedMacbethSamples(:,i) = cell2mat(cellfun(@nanmean, samples, 'UniformOutput',false)');
end


%% Load and crop a real image

fName = fullfile(parentPath,'Images','GOPR0049.ppm');
img = imread(fName);

img = double(img);
img = img/max(img(:));

x = 1422;
y = 1906;
delta = 10;

wp = img(y-delta:y+delta,x-delta:x+delta,:);
wp = reshape(wp,[(2*delta+1)^2, 3]);
wp = mean(wp);

imgVec = reshape(img,[size(img,1)*size(img,2),3]);
imgCorrVec = imgVec*diag(1./wp(:));
imgCorr = reshape(imgCorrVec,[size(img,1), size(img,2),3]);
imgCorr = imgCorr/max(imgCorr(:));

realCrop = imgCorr(967:2069,1299:2962,:);
realCrop = imresize(realCrop,[size(simCrop,1), size(simCrop,2)]);

realCrop = realCrop/max(realCrop(:));

figure; imshow(realCrop.^(1/2.2));

% Extract Macbeth patch intensities.
cp = [1350 2014;2894 2028;2900 997;1344 1010];

sensor = sensorCreate('monochrome');
sensor = sensorSet(sensor,'size',[size(imgCorr,1) size(imgCorr,2)]);
capturedMacbethSamples = zeros(24,3);

for i=1:3
   sensor = sensorSet(sensor,'volts',imgCorr(:,:,i));
   ieAddObject(sensor);
   sensorWindow();
   
   [samples, ~, ~, cp] = chartSelect(sensor,1,1,4,6,cp);
   capturedMacbethSamples(:,i) = cell2mat(cellfun(@nanmean, samples, 'UniformOutput',false)');
end


%% Create the plot

fig = figure;
hold on; grid on; box on;
plot(data.wave,data.cameraIll(:,1),'r','lineWidth',2);
plot(data.wave,data.cameraIll(:,2),'g','lineWidth',2);
plot(data.wave,data.cameraIll(:,3),'b','lineWidth',2);

plot(data.wave,data.ledCameraApprox(:,1),'r--','lineWidth',3);
plot(data.wave,data.ledCameraApprox(:,2),'g--','lineWidth',3);
plot(data.wave,data.ledCameraApprox(:,3),'b--','lineWidth',3);
xlabel('Wavelength, nm');
set(gca,'XTick',400:100:800);
set(gca,'FontSize',10);
set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',[1 1 9 3]);

axes('position',[0.65 0.6 0.3 0.3]);
annotation(fig,'textbox',...
    [0.75 0.55 0.108928571428571 0.05],'String',{'Captured'},'FontSize',10,'linestyle','none');
imshow(realCrop.^(1/2.2));

axes('position',[0.65 0.2 0.3 0.3]);
annotation(fig,'textbox',...
    [0.75 0.15 0.108928571428571 0.05],'String',{'Emulated'},'FontSize',10,'linestyle','none');
imshow(simCrop.^(1/2.2));

if isempty(destDir) == 0
    fName = fullfile(destDir,'RGBEmulationCurvesV2.eps');
    print('-depsc',fName);
end

% Scatter plot of emulated vs. captured Macbeth patch intensities.

figure;
hold on; grid on; box on;
plot(emulatedMacbethSamples(:,1),capturedMacbethSamples(:,1),'r.','markerSize',20);
plot(emulatedMacbethSamples(:,2),capturedMacbethSamples(:,2),'g.','markerSize',20);
plot(emulatedMacbethSamples(:,3),capturedMacbethSamples(:,3),'b.','markerSize',20);
xlabel('Emulated data');
ylabel('Captured data');






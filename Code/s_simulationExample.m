% This is a script illustrating the idea behind optimal illumination for
% pixel classification.
%
% In this script we use the ISET simulation environment to create a test
% chart composed of reflectances of different teeth samples. The samples
% were taken by measuring the reflectance spectra of a VITA Toothguide
% 3D-MASTER.
%
% We simulate data acquisition in four configurations:
% 1. A conventional RGB camera with a broadband 6500K illuminant. (To avoid
%    demosaicing issues we simulate each of the RGB channels separately).
% 2. A monochrome camera and 3 Bayer pattern-like illuminants (red, green
%    and blue).
% 3. A monochrome camera and 3 optimal illuminants estimated with the
%    unsupervised method.
% 4. A monochrome camera and 3 optimal illuminants estimated with the
%    supervised method.
%
% For each acquisition condition we run a SVM classifier with 10-fold cross
% validation on raw pixel data.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

ieInit;

[codePath, parentPath] = olRootPath();

wave = 400:4:1000;
nWaves = length(wave);

%% Initialization

% Load the spectra of the LEDs composing the LED cube.
fName = fullfile(parentPath,'Parameters','LEDCubeSpectra');
leds = ieReadSpectra(fName,wave);
leds = Energy2Quanta(wave,leds);
leds = leds/max(leds(:));

% Load the quantum efficiency of a silicone sensor. 
fName = fullfile(parentPath,'Parameters','qe');
qe = ieReadSpectra(fName,wave);

% Load the reflectances of teeth samples.
fName = fullfile(parentPath,'Parameters','teeth');
refl = ieReadSpectra(fName,wave);
reflCell = mat2cell(refl,nWaves,ones(size(refl,2),1));
reflCell = reshape(reflCell(1:28),[4,7]);

% Create a standard daylight-like, blackbody radiator.
ill = illuminantCreate('blackbody',wave,6500);
illQuanta = illuminantGet(ill,'photons');

% We assume that any light can only be generated with LEDCube, hence we
% need to find the best approximation of the broadband light
illApprox = leds*(leds\illQuanta);

figure;
hold on; grid on; box on;
plot(wave,[illApprox illQuanta]);
xlabel('Wavelength, nm');
ylabel('Illuminant');
title('Broadband illuminant approximation');
legend('LED approx.','original');

% Create an ISET scene template which contains a 4x7 chart of teeth
% reflectances.
sceneTemplate = sceneFromReflectanceCell(reflCell,10,wave);
sceneTemplate = sceneSet(sceneTemplate,'fov',10);
hfov = sceneGet(sceneTemplate,'fov horizontal');
vfov = sceneGet(sceneTemplate,'fov vertical');
sceneTemplate = sceneAdjustIlluminant(sceneTemplate,Quanta2Energy(wave,illApprox),0);
ieAddObject(sceneTemplate);
sceneWindow;

sensorTemplate = sensorCreate('monochrome');
sensorTemplate = sensorSet(sensorTemplate,'wave',wave);
sensorTemplate = sensorSet(sensorTemplate,'quantizationmethod','8bit');

%% 1. Conventional color camera under 6500K light

fName = fullfile(parentPath,'Parameters','GoProHero5.mat');
camera = ieReadColorFilter(wave,fName);

measVals = cell(1,28);
cameraGain = zeros(size(camera,2),1);

% Simulate data acquisition for each channel independently
for i=1:size(camera,2)
    
    oi = oiCreate;
    oi = oiCompute(oi,sceneTemplate);
    
    sensor = sensorSet(sensorTemplate,'filter spectra',camera(:,i));
    sensor = sensorSetSizeToFOV(sensor,[hfov vfov],sceneTemplate,oi);
    
    sensor = sensorCompute(sensor,oi);
    
    ieAddObject(sensor);
    sensorWindow;
    
    cameraGain(i) = sensorGainAndOffset(sceneTemplate,oi,sensor);
    
    sz = sensorGet(sensor,'size');
    cp = [1 sz(1); sz(2) sz(1); sz(2) 1; 1 1];
    
    tmp = chartSelect(sensor,1,1,4,7,cp);
    measVals = cellfun(@(x,y) cat(2,x,y),measVals,tmp,'UniformOutput',false);
end

% Assemble the data and run a classifier
nExamples = size(measVals{1},1);

data = cell2mat(measVals')*diag(cameraGain/max(cameraGain));
labels = repmat(1:28,[nExamples 1]);
labels = labels(:);

c = logspace(-2,2,10);
accy = zeros(length(c),1);
for b=1:length(c)
    template = templateSVM('Standardize',0,...
                           'KernelFunction','linear',...
                           'SaveSupportVectors','on',...
                           'BoxConstraint',c(b));

    model = fitcecoc(data,labels,'Learners',template,'Kfold',10);
    % cvModel = crossval(model,'Kfold',10);
    accy(b) = 1 - kfoldLoss(model);
    fprintf('RGB camera: c=%.4f, accuracy=%.2f\n',c(b),accy(b));
end
    
figure;
semilogx(c,accy);
grid on; box on;
ylim([0 1]);
xlabel('SVM c parameter');
ylabel('10-fold classification accy.');
title('Conventional + 6500K');


%% 2. Monochrome camera with Bayer pattern-like illuminants

% We generate illuminants that approximate the spectral responsivity
% functions of a Bayer pattern.
cameraApprox = leds*(leds\camera);

figure;
hold on; grid on; box on;
plot(wave,camera);
plot(wave,cameraApprox,'--');
xlabel('Wavelength, nm');
ylabel('Quanta au.');
title('Lights approximating Bayer filter transmissivities');
legend('Original','Approx.');


measVals = cell(1,28);
cameraGain = zeros(size(cameraApprox,2),1);

for i=1:size(cameraApprox,2)
    
    ill = cameraApprox(:,i)*10^16;
    
    scene = sceneAdjustIlluminant(sceneTemplate,Quanta2Energy(wave,ill),0);
    ieAddObject(scene);
    sceneWindow;
    
    oi = oiCreate;
    oi = oiCompute(oi,scene);

    sensor = sensorSet(sensorTemplate,'pixel spectral qe',ones(nWaves,1));
    sensor = sensorSetSizeToFOV(sensor,[hfov vfov],sceneTemplate,oi);
    
    sensor = sensorCompute(sensor,oi);
    
    cameraGain(i) = sensorGainAndOffset(sceneTemplate,oi,sensor);
    
    ieAddObject(sensor);
    sensorWindow;
    
    sz = sensorGet(sensor,'size');
    cp = [1 sz(1); sz(2) sz(1); sz(2) 1; 1 1];
    
    tmp = chartSelect(sensor,1,1,4,7,cp);
    measVals = cellfun(@(x,y) cat(2,x,y),measVals,tmp,'UniformOutput',false);
    
end

% Assemble the data and run a classifier
nExamples = size(measVals{1},1);
data = cell2mat(measVals')*diag(cameraGain/max(cameraGain));

labels = repmat(1:28,[nExamples 1]);
labels = labels(:);

c = logspace(-2,2,10);
accy = zeros(length(c),1);

for b=1:length(c)
    template = templateSVM('Standardize',0,...
                           'KernelFunction','linear',...
                           'SaveSupportVectors','on',...
                           'BoxConstraint',c(b));

    model = fitcecoc(data,labels,'Learners',template,'Kfold',10);
    accy(b) = 1 - kfoldLoss(model);
    fprintf('RGB-flash camera: c=%.4f, accuracy=%.2f\n',c(b),accy(b));
end
    
figure;
semilogx(c,accy);
grid on; box on;
ylim([0 1]);
xlabel('SVM c parameter');
ylabel('10-fold classification accy.');
title('RGB-flash like camera');


%% Optimal illumination estimated using the unsupervised approach

lambda = 0.01;

% We need to account for sensor qe when we compute the optimal illuminant
dta = diag(qe)*refl;
X = sparsePCA(cov(dta'), 3, lambda, 1e-6, leds, ...
    'maxIter',100, 'Verbose', 0);

% Sometimes X may contain very small negative numbers.
X = max(X,0);


figure;
hold on; grid on; box on;
plot(wave,X,'LineWidth',2);
xlabel('Wavelength, nm');
ylabel('Scaled quanta');
title('Optimal illuminants (unsupervised)');

measVals = cell(1,28);
cameraGain = zeros(size(X,2),1);

for i=1:size(X,2)
    
    ill = leds*(leds\X(:,i))*10^16;
    
    scene = sceneAdjustIlluminant(sceneTemplate,Quanta2Energy(wave,ill),0);
    ieAddObject(scene);
    sceneWindow;
    
    oi = oiCreate;
    oi = oiCompute(oi,scene);

    
    sensor = sensorSet(sensorTemplate,'pixel spectral qe',qe);
    sensor = sensorSetSizeToFOV(sensor,[hfov vfov],sceneTemplate,oi);
    
    sensor = sensorCompute(sensor,oi);
    
    cameraGain(i) = sensorGainAndOffset(sceneTemplate,oi,sensor);
    
    ieAddObject(sensor);
    sensorWindow;
    
    sz = sensorGet(sensor,'size');
    cp = [1 sz(1); sz(2) sz(1); sz(2) 1; 1 1];
    
    tmp = chartSelect(sensor,1,1,4,7,cp);
    measVals = cellfun(@(x,y) cat(2,x,y),measVals,tmp,'UniformOutput',false);
    
end


% Assemble the data and run a classifier
nExamples = size(measVals{1},1);
data = cell2mat(measVals')*diag(cameraGain/max(cameraGain));

labels = repmat(1:28,[nExamples 1]);
labels = labels(:);

c = logspace(-2,2,10);
accy = zeros(length(c),1);
for b=1:length(c)
    template = templateSVM('Standardize',0,...
                           'KernelFunction','linear',...
                           'SaveSupportVectors','on',...
                           'BoxConstraint',c(b));

    model = fitcecoc(data,labels,'Learners',template,'Kfold',10);
    accy(b) = 1 - kfoldLoss(model);
    fprintf('Unsupervised illuminants: c=%.4f, accuracy=%.2f\n',c(b),accy(b));
end
    
figure;
semilogx(c,accy);
grid on; box on;
ylim([0 1]);
xlabel('SVM c parameter');
ylabel('10-fold classification accy.');
title('Unsupervised illuminants');

%% Optimal illumination estimated using the supervised approach


dta = diag(qe)*refl;
dta = dta/max(dta(:));
[ phi, predLabels, W, b, hist ] = multiclassSVM( dta, 1:size(dta,2) , leds,...
            'multiClassMaxIter', 10, 'multiclassC', 10);

figure;
hold on; grid on; box on;
plot(wave,phi);
xlabel('Wavelength, nm');
ylabel('Scaled quanta');
title('Optimal illuminants (supervised)');

measVals = cell(1,28);
cameraGain = zeros(size(phi,2),1);

for i=1:size(X,2)
    
    ill = leds*(leds\phi(:,i))*10^16;
    
    scene = sceneAdjustIlluminant(sceneTemplate,Quanta2Energy(wave,ill),0);
    ieAddObject(scene);
    sceneWindow;
    
    oi = oiCreate;
    oi = oiCompute(oi,scene);

    
    sensor = sensorSet(sensorTemplate,'pixel spectral qe',qe);
    sensor = sensorSetSizeToFOV(sensor,[hfov vfov],sceneTemplate,oi);
    
    sensor = sensorCompute(sensor,oi);
    
    cameraGain(i) = sensorGainAndOffset(sceneTemplate,oi,sensor);
    
    ieAddObject(sensor);
    sensorWindow;
    
    sz = sensorGet(sensor,'size');
    cp = [1 sz(1); sz(2) sz(1); sz(2) 1; 1 1];
    
    tmp = chartSelect(sensor,1,1,4,7,cp);
    measVals = cellfun(@(x,y) cat(2,x,y),measVals,tmp,'UniformOutput',false);
    
end


% Assemble the data and run a classifier
nExamples = size(measVals{1},1);
data = cell2mat(measVals')*diag(cameraGain/max(cameraGain));

labels = repmat(1:28,[nExamples 1]);
labels = labels(:);

c = logspace(-2,2,10);
accy = zeros(length(c),1);
for b=1:length(c)
    template = templateSVM('Standardize',0,...
                           'KernelFunction','linear',...
                           'SaveSupportVectors','on',...
                           'BoxConstraint',c(b));

    model = fitcecoc(data,labels,'Learners',template,'Kfold',10);
    accy(b) = 1 - kfoldLoss(model);
    fprintf('Supervised illuminants: c=%.4f, accuracy=%.2f\n',c(b),accy(b));
end
    
figure;
semilogx(c,accy);
ylim([0 1]);
xlabel('SVM c parameter');
ylabel('10-fold classification accy.');
title('Supervised illuminants');


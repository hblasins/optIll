% This is a script illustrating how the optimal selection algorithms can
% help understand the tradeoff between the number of acquisitions and pixel
% classification accuracy.
%
% In this script we use the ISET simulation environment to create a test
% chart composed of reflectances of different teeth samples. The samples
% were taken by measuring the reflectance spectra of a VITA Toothguide
% 3D-MASTER.
%
% We simulate data acquisition in three configurations:
% 1. A conventional RGB camera with a broadband 6500K illuminant. (To avoid
%    demosaicing issues we simulate each of the RGB channels separately). 
% 2. A monochrome camera capturing data under all LEDCube illuminants
%    individually switched on (11 channels in total).
% 3. A monochrome camera and 3 optimal illuminants estimated with the
%    unsupervised method.
%
% For each acquisition condition we run a SVM classifier and evaluate its
% performance with a 30% held-out set (this is faster than 10-fold method).
%
% Copyright, Henryk Blasinski 2017.

close all;
clear variables;
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
nLeds = size(leds,2);

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


sceneTemplate = sceneFromReflectanceCell(reflCell,10,wave);
sceneTemplate = sceneSet(sceneTemplate,'fov',10);
hfov = sceneGet(sceneTemplate,'fov horizontal');
vfov = sceneGet(sceneTemplate,'fov vertical');
ieAddObject(sceneTemplate);
sceneWindow;

sensorTemplate = sensorCreate('monochrome');
sensorTemplate = sensorSet(sensorTemplate,'wave',wave);
sensorTemplate = sensorSet(sensorTemplate,'pixel spectral qe',qe');
sensorTemplate = sensorSet(sensorTemplate,'quantizationmethod','8bit');

%% 1. Conventional RGB camera and broadband 6500K illuminant

fName = fullfile(parentPath,'Parameters','GoProHero5.mat');
camera = ieReadColorFilter(wave,fName);

% Simulate data acquisition for each channel independently
measVals = cell(1,28);
cameraGain = zeros(size(camera,2),1);
for i=1:size(camera,2)
    
    scene = sceneAdjustIlluminant(sceneTemplate,Quanta2Energy(wave,illApprox),0);
    
    oi = oiCreate;
    oi = oiCompute(oi,scene);
    
    sensor = sensorSet(sensorTemplate,'filter spectra',camera(:,i));
    sensor = sensorSetSizeToFOV(sensor,[hfov vfov],scene,oi);
    sensor = sensorCompute(sensor,oi);
    
    ieAddObject(sensor);
    sensorWindow;
    
    cameraGain(i) = sensorGainAndOffset(scene,oi,sensor);
    
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

    model = fitcecoc(data,labels,'Learners',template,'Holdout',0.3);
    accy(b) = 1 - kfoldLoss(model);
    fprintf('RGB camera + blackbody 6500K: c=%.4f, accuracy=%.2f\n',c(b),accy(b));
end
    
figure;
semilogx(c,accy);
grid on; box on;
ylim([0 1]);
xlabel('SVM c parameter');
ylabel('10-fold classification accy.');
title('Conventional + 6500K');

conventionalAccuracy = max(accy);

%% 2. Monochrome camera and multispectral flash capture (11 LEDs)

measVals = cell(1,28);
cameraGain = zeros(nLeds,1);
for i=1:nLeds
    
    ill = leds(:,i)*10^16;
    
    scene = sceneAdjustIlluminant(sceneTemplate,Quanta2Energy(wave,ill),0);
    scene = sceneSet(scene,'name',sprintf('Teeth + LED %i',i));
    ieAddObject(scene);
    sceneWindow;
    
    oi = oiCreate;
    oi = oiCompute(oi,scene);
  
    sensor = sensorSetSizeToFOV(sensorTemplate,[hfov vfov],sceneTemplate,oi);    
    sensor = sensorCompute(sensor,oi);
    
    cameraGain(i) = sensorGainAndOffset(sceneTemplate,oi,sensor);
    
    ieAddObject(sensor);
    sensorWindow;
    
    sz = sensorGet(sensor,'size');
    cp = [1 sz(1); sz(2) sz(1); sz(2) 1; 1 1];
    
    tmp = chartSelect(sensor,1,1,4,7,cp);
    measVals = cellfun(@(x,y) cat(2,x,y),measVals,tmp,'UniformOutput',false);
    
end

% Extract the data and run classification algorithm
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

    model = fitcecoc(data,labels,'Learners',template,'Holdout',0.3);
    accy(b) = 1 - kfoldLoss(model);
    fprintf('Multispectral flash: c=%.4f, accuracy=%.2f\n',c(b),accy(b));
end

multispectralAccuracy = max(accy);

%% 3. Optimal illumination estimated using the unsupervised approach
%  We vary the number of optimal lights between 1 and 11 to evaluate the
%  performance gains.

Kmax = 5;

optimalIlluminants = cell(1,Kmax);
optimalAccuracy = zeros(1,Kmax);

for kk=1:Kmax

    % The spectral covariance matrix needs to account for sensor quantum
    % efficiency
    dta = diag(qe)*refl;
    [X, wghts, S, hist] = sparsePCA(cov(dta'), kk, 0.005, 1e-6, leds, ...
        'maxIter',100, ...
        'Verbose', 0);
    X = max(X,0);

    optimalIlluminants{kk} = X;

    figure;
    hold on; grid on; box on;
    plot(wave,X,'LineWidth',2);
    xlabel('Wavelength, nm');
    ylabel('Optimal illuminant spectra');
    title(sprintf('Number of optimal lights: %i ',kk));
    drawnow;
    
    % Simulate the acquisition
    measVals = cell(1,28);
    cameraGain = zeros(kk,1);
    for i=1:kk
        
        ill = leds*(leds\X(:,i))*10^16;
        
        scene = sceneAdjustIlluminant(sceneTemplate,Quanta2Energy(wave,ill),0);
        ieAddObject(scene);
        sceneWindow;
        
        oi = oiCreate;
        oi = oiCompute(oi,scene);
        
        
        sensor = sensorSetSizeToFOV(sensorTemplate,[hfov vfov],scene,oi);
        sensor = sensorCompute(sensor,oi);
        ieAddObject(sensor);
        sensorWindow;
        
        cameraGain(i) = sensorGainAndOffset(scene,oi,sensor);
        
        sz = sensorGet(sensor,'size');
        cp = [1 sz(1); sz(2) sz(1); sz(2) 1; 1 1];
        
        tmp = chartSelect(sensor,1,1,4,7,cp);
        measVals = cellfun(@(x,y) cat(2,x,y),measVals,tmp,'UniformOutput',false);
        
    end


    % Assemble the data and run classification algorithm
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
        
        model = fitcecoc(data,labels,'Learners',template,'Holdout',0.3);
        accy(b) = 1 - kfoldLoss(model);
        fprintf('Optimal flash (Unsup.), %i flash(es): c=%.4f, accuracy=%.2f\n',kk,c(b),accy(b));
    end
    
    optimalAccuracy(kk) = max(accy);

end

%% Compare approaches

figure;
hold on; grid on; box on;
plot(conventionalAccuracy*ones(Kmax,1),':r','LineWidth',2);
plot(multispectralAccuracy*ones(Kmax,1),'--g','LineWidth',2);
plot(optimalAccuracy,'LineWidth',2);
xlabel('Number of illuminants');
ylabel('Accuracy');
ylim([0.5,1]);
set(gca,'XTick',1:Kmax);
legend('RGB + 6500K','Multispctral','Optimal (Unsup.)','Location','SouthEast');



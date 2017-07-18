% This script builds a synthetic example of a scene where the feature of
% interest (an ellipse) is invisible under broadband illuminant. However
% camera is able to detect it, once an optimal illuminant is estimated and
% used.
%
% This toy example illustrates the key idea and the usage of the
% algorithms.
%
% Copyright, Henryk Blasinski 2017
close all;
clear all;
clc;

ieInit;

[codePath, parentPath] = olRootPath();

wave = 400:10:700;
nWaves = length(wave);

%% Initialize

% Load a camera and some example reflectance spectra
fName = fullfile(parentPath,'Parameters','GoProHero5');
camera = ieReadColorFilter(wave,fName);

% Load Macbeth reflectance data
fName = fullfile(parentPath,'Parameters','macbethChart');
refl = ieReadSpectra(fName,wave);

% Pick a reflectance chip
chipID = 7;

%% Build two different reflectances that are metamers to the camera
%  Search for the smoothest spectrum in that gives the same camera response
%  as the reference spectrum.
R = [eye(nWaves - 1), zeros(nWaves-1,1)] - [zeros(nWaves-1,1), eye(nWaves - 1)];

cvx_begin
    variable spec(nWaves,1)
    minimize norm(R*spec)
    subject to
        1 >= spec >= 0
        camera'*spec == camera'*refl(:,chipID)
cvx_end

figure;
hold on; grid on; box on;
plot(wave,[spec, refl(:,chipID)]);
xlabel('Wavelength, nm');
ylabel('Reflectance');

% Create a scene out of the two spectra (an ellipse)

hh = 100;
ww = 100;
mask = zeros(hh,ww);
for xx=1:ww
    for yy=1:hh
        if (xx-50)^2 + ((yy-50)/2)^2 <= 300
            mask(yy,xx) = 1;
        end
    end
end

figure;
imshow(mask,'Border','tight','InitialMagnification',500);
colormap gray; axis image;

bbandLight= illuminantCreate('blackbody',wave,6500);
bbandLight = illuminantGet(bbandLight,'photons');
bbandLight = bbandLight/max(bbandLight(:));



% Create a scene
reflData = zeros(nWaves,hh*ww);
reflData(:,find(mask(:) > 0)) = repmat(diag(bbandLight)*refl(:,chipID),[1 sum(mask(:))]);
reflData(:,find(mask(:)==0)) = repmat(diag(bbandLight)*spec,[1 sum(mask(:)==0)]);
reflData = reflData*10^15;

reflCell = mat2cell(reflData,nWaves,ones(1,hh*ww));
reflCell = reshape(reflCell,[hh ww]);

scene = sceneFromReflectanceCell(reflCell,1,wave);
scene = sceneSet(scene,'name','Broadband illuminant');
ieAddObject(scene);

% Create an oi
oi = oiCreate;
oi = oiCompute(oi,scene);
oi = oiSet(oi,'name','Broadband illuminant');
ieAddObject(oi);

% Create a sensor
sensor = sensorCreate('bayer (rggb)');
sensor = sensorSet(sensor,'size',[hh ww]);
sensor = sensorSet(sensor,'filter transmissivities',camera);
sensor = sensorSet(sensor,'noise flag',0);
sensor = sensorSetSizeToFOV(sensor,sceneGet(scene,'fov'),scene,oi);
sensor = sensorSet(sensor,'name','Broadband illuminant');
sensor = sensorCompute(sensor,oi);
ieAddObject(sensor);

% Create the digital image
ip = ipCreate;
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name','Broadband illuminant');
ieAddObject(ip);

sceneWindow;
oiWindow;
sensorWindow;
ipWindow;


%% Compute an optimal illuminant (unsupervised approach)

fName = fullfile(parentPath,'Parameters','LEDCubeSpectra');
leds = ieReadSpectra(fName,wave);
leds = Energy2Quanta(wave,leds);
leds = leds/max(leds(:));

dta = [spec, refl(:,chipID)];
X = sparsePCA(cov(dta'), 1, 0.001, 1e-6, leds, 'maxIter',100, 'Verbose', 0);

figure;
hold on; grid on; box on;
plot(wave,X/max(X(:)),'b','lineWidth',2);
xlabel('Wavelength, nm');
title('Optimal illuminant');
ylim([-0.05 1.05]);

% Illuminate the scene with the optimal illuminant
sceneOptimal = sceneAdjustIlluminant(scene,Quanta2Energy(wave,X*10^15),0);
sceneOptimal = sceneSet(sceneOptimal,'name','Optimal illuminant');
ieAddObject(sceneOptimal);

% Optical image
oiOptimal = oiCompute(oi,sceneOptimal);
oiOptimal = oiSet(oiOptimal,'name','Optimal illuminant');
ieAddObject(oiOptimal);

% Sensor
sensorOptimal = sensorCompute(sensor,oiOptimal);
sensorOptimal = sensorSet(sensorOptimal,'name','Optimal illuminant');
ieAddObject(sensorOptimal);

% ISP
ipOptimal = ipCompute(ip,sensorOptimal);
ipOptimal = ipSet(ipOptimal,'name','Optimal illuminant');
ieAddObject(ipOptimal);

sceneWindow;
oiWindow;
sensorWindow;
ipWindow;




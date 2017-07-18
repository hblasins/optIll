function [ scene ] = sceneFromReflectanceCell(reflectanceCell,pSize,wave)
%
% Create a target with reflectances specified in the reflectanceCell and
% preserving their spatial arrangement.
%
% Copyright Henryk Blasinski, 2014

if ieNotDefined('pSize'),    pSize = 32; end

% Default scene
scene = sceneCreate;
if ieNotDefined('wave'), wave = sceneGet(scene,'wave');
else                     scene = sceneSet(scene,'wave',wave);
end
nWave = length(wave);
defaultLuminance = 100;  % cd/m2

% Spatial arrangement
r = size(reflectanceCell,1); c = size(reflectanceCell,2);
rcSize = [r,c];

% Convert the scene reflectances into photons assuming an equal photon
% illuminant.
ee         = ones(nWave,1);           % Equal photon vector

% Put these into the scene data structure.  These are in photon units, but
% they are not scaled to reasonable photon values.
sData = zeros(rcSize(1),rcSize(2),nWave);
for rr=1:rcSize(1)
    for cc=1:rcSize(2)
        refl = reflectanceCell{rr,cc};
        if isempty(refl), refl = zeros(nWave,1); end;
        sData(rr,cc,:) = diag(ee)*squeeze(refl);
    end
end

% Build up the size of the image regions - still reflectances
sData = imageIncreaseImageRGBSize(sData,pSize);

% Add data to scene, using equal energy illuminant
scene = sceneSet(scene,'photons',sData);
scene = sceneSet(scene,'illuminantPhotons',ee);
scene = sceneSet(scene,'illuminantComment','Equal photon');
scene = sceneSet(scene,'name','Reflectance Chart (EE)');
% vcAddAndSelectObject(scene); sceneWindow;

% Adjust the illuminance to a default level in cd/m2
% scene = sceneAdjustLuminance(scene,defaultLuminance);
% vcAddAndSelectObject(scene); sceneWindow;

return




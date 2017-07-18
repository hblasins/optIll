% This script controlls the LEDcube illumination and the PointGrey camera
% to capture images of a scene using the RGB camera emulation approach.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear variables;
clc;

ieInit;

[codePath, parentPath] = olRootPath();

% Assign a name to the scene.
target = 'Lemons';

wave = 400:4:800;
nWaves = length(wave);

%% Initialization

fName = fullfile(parentPath,'Parameters','LEDCubeSpectra');
leds = ieReadSpectra(fName,wave);
leds = Energy2Quanta(wave,leds);
leds = leds/max(leds(:));
nLEDs = size(leds,2);


fName = fullfile(parentPath,'Parameters','qe');
qe = ieReadSpectra(fName,wave);

cameras = {'AptinaMT9M031','AptinaMT9M131',...
    'Canon1DMarkIII','Canon5DMarkII','Canon20D','Canon40D','Canon50D','Canon60D','Canon300D','Canon500D','Canon600D'...
    'HasselbladH2',...
    'NikonD1','NikonD3','NikonD3X','NikonD40','NikonD50','NikonD70','NikonD80','NikonD90','NikonD100','NikonD200',...
    'NikonD200IR','NikonD300s','NikonD700','NikonD5100',...
    'NokiaN900',...
    'OlympusE-PL2',...
    'PentaxK-5','PentaxQ',...
    'PhaseOne',...
    'PointGreyGrasshopper50S5C','PointGreyGrasshopper214S5C',...
    'SONYNEX-5N'};

illTemps = [10000 6500 4000 2000];

%% Main acquisition loop

for j=1:length(illTemps)
    
    % Generate the desited illuminant spectrum.
    ill = illuminantCreate('blackbody',wave,illTemps(j));
    ill = illuminantGet(ill,'photons');
    ill = ill/max(ill);
    
    for n=1:length(cameras)
        
        
        % Load camera responsivity functions
        referenceCamera = cameras{n};
        fName = fullfile(parentPath,'Parameters','Cameras',referenceCamera);
        camera = ieReadColorFilter(wave,fName);
        camera(isnan(camera)) = 0;
        
       
        % Compute the illuminant that emulates the combination of the
        % desired illuminant-camera pair.
        cvx_begin
            variable ledCameraApproxWghts(nLEDs,3)
            minimize sum(norms(diag(ill)*camera - diag(qe)*leds*ledCameraApproxWghts,2,1))
            subject to
                ledCameraApproxWghts >= 0
                % We omit wideband LEDs from the optimization to improve
                % color rendering. The norm objective in wavelength space
                % is not the same as color rendering.
                ledCameraApproxWghts([1, 3, 5],:) == 0

        cvx_end
        
        cameraIll = diag(ill)*camera;
        ledCameraApprox = diag(qe)*leds*ledCameraApproxWghts;
        
        
        figure;
        hold on; grid on; box on;
        plot(wave,cameraIll,'LineWidth',2);
        plot(wave,ledCameraApprox,'--','LineWidth',2);
        xlabel('Wavelength, nm');
        ylabel('Scaled quanta');
        title(sprintf('%s + %iK emulation',referenceCamera,illTemps(j)));
        
       
        %% Capture a set of images as if using a monochrome camera
        %  and the lights that emulate the defined RGB camera.
        
        Img = zeros(1024,1280,3);
        ImgRAW = zeros(1024,1280,3);
        sh = zeros(3,1);
        gn = zeros(3,1);
        for i=1:3
            [ImgRAW(:,:,i), sh(i), gn(i)] = getFleaLEDCube( -1, -1, 0, 10, 'Tmp', 10, 0, ledCameraApproxWghts(:,i) );
            
            % We store data in the linear space (i.e. we account for gain
            % and exposure duraion).
            linGain = sh(i)*10^(gn(i)/20);
            Img(:,:,i) = ImgRAW(:,:,i)/linGain;
        end
              
        % But we normalize image intensities.
        Img = Img/max(Img(:));
        
        fName = fullfile(parentPath,'Images',target,'Conventional',sprintf('%s_%iK.mat',referenceCamera,illTemps(j)));
        dirName = fileparts(fName);
        if ~exist(dirName,'dir'),
            mkdir(dirName);
        end
        
        save(fName,'Img','ledCameraApprox','cameraIll','wave');
        
    end
end

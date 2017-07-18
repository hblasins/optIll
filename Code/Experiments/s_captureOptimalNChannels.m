% This script controlls the LEDcube illumination and the PointGrey camera
% to capture images of a scene using the illuminants computed via the
% supervised and unsupervised approaches. The selection algorithms are
% cross-validated for their respective tuning parameters as well as the
% number of illuminants.
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

%%

fName = fullfile(parentPath,'Parameters','LEDCubeSpectra');
leds = ieReadSpectra(fName,wave);
leds = Energy2Quanta(wave,leds);
leds = leds/max(leds(:));
nLEDs = size(leds,2);

fName = fullfile(parentPath,'Parameters','qe');
qe = ieReadSpectra(fName,wave);

% For every scene you need to provide a set of reflectance spectra and the
% corresponding class labels.
switch target
        
    case 'GreenApples'
        
        fName = fullfile(olRootPath,'Data','realGreenApple');
        refl1 = ieReadSpectra(fName,wave);
        labels1 = ones(1,size(refl1,2));
        
        fName = fullfile(olRootPath,'Data','fakeGreenApple');
        refl2 = ieReadSpectra(fName,wave);
        labels2 = 2*ones(1,size(refl2,2));
        
        refl = [refl1, refl2];
        labels = [labels1(:); labels2(:)];
        
    case 'Lemons'
        
        fName = fullfile(olRootPath,'Data','realLemon');
        refl1 = ieReadSpectra(fName,wave);
        labels1 = ones(1,size(refl1,2));
        
        fName = fullfile(olRootPath,'Data','fakeLemon');
        refl2 = ieReadSpectra(fName,wave);
        labels2 = 2*ones(1,size(refl2,2));
        
        refl = [refl1, refl2];
        labels = [labels1(:); labels2(:)];
        
    case 'RedAndGreenApples'
        
        fName = fullfile(olRootPath,'Data','realGreenApple');
        refl1 = ieReadSpectra(fName,wave);
        labels1 = ones(1,size(refl1,2));
        
        fName = fullfile(olRootPath,'Data','fakeGreenApple');
        refl2 = ieReadSpectra(fName,wave);
        labels2 = 2*ones(1,size(refl2,2));
        
        fName = fullfile(olRootPath,'Data','realRedApple');
        refl3 = ieReadSpectra(fName,wave);
        labels3 = 3*ones(1,size(refl3,2));
        
        fName = fullfile(olRootPath,'Data','fakeRedApple');
        refl4 = ieReadSpectra(fName,wave);
        labels4 = 4*ones(1,size(refl4,2));
        
        refl = [refl1, refl2, refl3, refl4];
        labels = [labels1(:); labels2(:); labels3(:); labels4(:)];
        
    case 'YellowAndGreenPears'
        
        fName = fullfile(olRootPath,'Data','realGreenPear');
        refl1 = ieReadSpectra(fName,wave);
        labels1 = ones(1,size(refl1,2));
        
        fName = fullfile(olRootPath,'Data','fakeGreenPear');
        refl2 = ieReadSpectra(fName,wave);
        labels2 = 2*ones(1,size(refl2,2));
        
        fName = fullfile(olRootPath,'Data','realYellowPear');
        refl3 = ieReadSpectra(fName,wave);
        labels3 = 3*ones(1,size(refl3,2));
        
        fName = fullfile(olRootPath,'Data','fakeYellowPear');
        refl4 = ieReadSpectra(fName,wave);
        labels4 = 4*ones(1,size(refl4,2));
        
        refl = [refl1, refl2, refl3, refl4];
        labels = [labels1(:); labels2(:); labels3(:); labels4(:)];
        
end


%% Unsupervised
%  Generate illuminants using the unsupervised approach. Cross-validate for
%  unsupervised sparsity parameters and the number of illuminants.
lambdas = logspace(-3,0,10);

for chnls=1:10
    for lambda=lambdas
        
        % First compute the optimal illuminant
        dta = diag(qe)*refl;
        [X, wghts, S, hist] = sparsePCA(cov(dta'), chnls, lambda, 1e-6, leds, 'maxIter',100, 'Verbose', 0);
        X = max(X,0);
        
        % We can only generate nonnegtive mixtures of LEDCube leds.
        cvx_begin
            variable wghts(nLEDs,chnls)
            minimize sum(norms(X - leds*wghts,2,1))
            subject to
                wghts >= 0
        cvx_end
        
        %{
        figure;
        hold on;
        plot(wave,X);
        plot(wave,leds*wghts,'--');
        xlabel('Wavelength, nm');
        title('Unsupervised illuminants');
        legend('Desired','Realizable');
        %}
        
        % Capture each channel independently
        Img = zeros(1024,1280,chnls);
        sh = zeros(chnls,1);
        gn = zeros(chnls,1);
        for i=1:chnls
            [Img(:,:,i), sh(i), gn(i)] = getFleaLEDCube( -1, 0, 0, 10, 'Tmp', 10, 0, wghts(:,i) );
            
            linGain = sh(i)*10^(gn(i)/20);
            Img(:,:,i) = Img(:,:,i)/linGain;
        end
        
        % Save data
        Img = Img/max(Img(:));
        
        fName = fullfile(parentPath,'Images',target,'Optimal',sprintf('%s_%f_%i.mat','Unsup',lambda,chnls));
        dirName = fileparts(fName);
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
        
        save(fName,'Img','X');
    end
end

%% Supervised
%  Generate illuminants using the supervised approach. Cross-validate for
%  multiclass SVM penalty and the number of illuminants.
lambdas = logspace(0,5,10);

for chnls=1:10
    for lambda=lambdas
        
        dta = diag(qe)*refl;
        [ phi, predLabels, W, b, hist ] = multiclassSVM( dta, labels' , leds,...
            'multiClassMaxIter', 10, 'multiclassC', lambda,...
            'nChannels', chnls, 'alpha', 1e-6);
        
        X = max(phi,0);
        
        % We can only generate nonnegtive mixtures of LEDCube leds.
        cvx_begin
            variable wghts(nLEDs,chnls)
            minimize sum(norms(X - leds*wghts,2,1))
            subject to
                wghts >= 0
        cvx_end
        
         %{
        figure;
        hold on;
        plot(wave,X);
        plot(wave,leds*wghts,'--');
        xlabel('Wavelength, nm');
        title('Supervised illuminants');
        legend('Desired','Realizable');
        %}
        
        Img = zeros(1024,1280,chnls);
        sh = zeros(chnls,1);
        gn = zeros(chnls,1);
        for i=1:chnls
            [Img(:,:,i), sh(i), gn(i)] = getFleaLEDCube( -1, 0, 0, 10, 'Tmp', 10, 0, wghts(:,i) );
            
            linGain = sh(i)*10^(gn(i)/20);
            Img(:,:,i) = Img(:,:,i)/linGain;
            
        end
        
        % Save data
        Img = Img/max(Img(:));
        
        fName = fullfile(parentPath,'Images',target,'Optimal',sprintf('%s_%f_%i.mat','Sup',lambda,chnls));
        dirName = fileparts(fName);
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
        
        save(fName,'Img','X');
    end
end
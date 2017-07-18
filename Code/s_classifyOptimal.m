% This is the main classification script for cameras using the optimized 
% illuminants. We evaluate illuminants generated with different tuning parameter 
% settings, as well as different number of optimal illuminants. Different
% classifiers are also cross-validated for tuning parameters. 
%
% Copyright, Henryk Blasinski 2017
close all;
clear all;
clc;

[codePath, parentPath] = olRootPath();
imagePath = fullfile(parentPath,'Images');
resultPath = fullfile(parentPath,'ResultsTest');

validTargets = {'GreenApples','YellowAndGreenPears','Lemons',...
    'RedAndGreenApples'};

target = 'Lemons';


%% Optimal cameras

fNames = fullfile(imagePath,target,'Optimal','*.mat');
fileNames = dir(fNames);

switch target
    
    case 'YellowAndGreenPears'
        x = [373;868;295;733];
        y = [241;221;593;509];
        delta = 50;
    
    case 'Lemons'
        x = [364; 890];
        y = [457; 454];
        delta = 100;
    
    case 'GreenApples'
        x = [300;910];
        y = [325;329];
        delta = 100;
        
    case 'RedAndGreenApples'
        x = [339;816;260;818];
        y = [257;263;531;552];
        delta = 50;
end

    
labels = zeros(length(x)*(2*delta+1)^2,1);   
for i=1:length(x)
    labels((i-1)*(2*delta+1)^2+1:i*(2*delta+1)^2) = i;
end

%% Main classification loop.
%  Classify for every conventional camera type, number of illuminants, cross
%  validate for different settings of validation parameters.

rng(123); % Reproducibility
partition = cvpartition(labels,'HoldOut',0.3);


classifierNames = {'SVM','KNN','DA','Tree','NB'};
regularizers = {logspace(-3,3,10),...
                round(linspace(1,50,10)),...
                linspace(0,1,10),...
                round(linspace(1,50,10)),...
                1};  

% delete(gcp('nocreate'))
% cluster = parcluster('local');
% cluster.NumWorkers = 32;
% pool = parpool(cluster,31);
            
for f=1:length(fileNames);
    
    fName = fullfile(imagePath,target,'Optimal',fileNames(f).name);
    [~, condition] = fileparts(fName);
    loadedData = load(fName);
    
    %{
    if f==1
        figure;
        imshow(loadedData.Img); 
        for r=1:length(x)
            rectangle('Position',[x(r)-delta y(r)-delta 2*delta+1 2*delta+1],'edgeColor','red');
        end
    end
    %}
    
    % Load and assemble the data.
    nChannels = size(loadedData.Img,3);
    data = zeros(length(x)*(2*delta+1)^2,nChannels);    
    for i=1:length(x)     
        subImg = loadedData.Img(y(i)-delta:y(i)+delta,x(i)-delta:x(i)+delta,:);
        subData = reshape(subImg,[(2*delta+1)^2 nChannels]);
        data((i-1)*(2*delta+1)^2+1:i*(2*delta+1)^2,:) = subData;        
    end    
    
    for c=1:length(classifierNames)
    
        fName = fullfile(resultPath,target,'Optimal',sprintf('%s_%s.mat',condition,classifierNames{c}));
        if exist(fName,'file')
            fprintf('File %s exists, skipping.\n',fName);
            continue;
        end
        
        reg = regularizers{c};
        
        accy = zeros(length(reg),1);
        predicted = cell(length(reg),1);
        true = cell(length(reg),1);
        model = cell(length(reg),1);
        
        for b=1:length(reg)
            
            switch classifierNames{c}
                case 'SVM'
                    template = templateSVM('Standardize',1,...
                        'KernelFunction','linear',...
                        'SaveSupportVectors','on',...
                        'BoxConstraint',reg(b));
                case 'KNN'
                    template = templateKNN('BreakTies','random',...
                                           'distance','euclidean',...
                                           'NumNeighbors',reg(b),...
                                           'Standardize',1);
                case 'DA'
                    template = templateDiscriminant('Gamma',reg(b));
                case 'Tree'
                    template = templateTree('MaxNumSplits',reg(b),...
                                           'Prune','on');
                case 'NB'
                    template = templateNaiveBayes();               
            end
            
            
            model{b} = fitcecoc(data,labels,'Learners',template,'CVPartition',partition);
            predLabels = kfoldPredict(model{b});
            
            predicted{b} = predLabels(test(partition));
            true{b} = labels(test(partition));
            
            accy(b) = 1 - kfoldLoss(model{b});
            fprintf('Camera: %s, Classifier: %s, reg.=%.4f, acc.=%.2f\n',condition,classifierNames{c},reg(b),accy(b));
        end
        
        
        if ~exist(fileparts(fName),'dir'),
            mkdir(fileparts(fName));
        end
        
        
        % parforSave(fName,accy,reg,predicted,true,model); 
        save(fName,'accy','reg','predicted','true','model');
    
    end
    
end

% delete(pool);



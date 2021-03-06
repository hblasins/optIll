% This script analyzes the classification results generated by s_classify.m
% script and investigates the relationship between the accuracy and the
% number of optimal illuminants. You can use this script to reproduce Fig.
% 4 from the manuscript.
%
% Note: this script takes a while to execute.
%
% Copyright, Henryk Blasinski 2017
close all;
clear all;
clc;

[codePath, parentPath] = olRootPath();

resultDir = fullfile(parentPath,'ResultsFull');
% resultDir = fullfile('/','Volumes','MyPassport','OptimalLight','ResultsSimplified');
% destPath = fullfile(parentPath,'TestFigures');
destPath = [];

target = 'RedAndGreenApples';

optimalCameras = {'Unsup','Sup'};
classifierNames = {'SVM','KNN','DA','Tree','NB'};

maxChnls = 10;

%% Read in classification results

accy = zeros(length(optimalCameras),length(classifierNames),maxChnls);
curves = cell(length(optimalCameras),length(classifierNames),maxChnls);

for o=1:length(optimalCameras)
    for c=1:length(classifierNames)
        
        for chnls=1:maxChnls
            
            inputDir = fullfile(resultDir,target,'Optimal',sprintf('%s_*_%i_%s.mat',optimalCameras{o},chnls,classifierNames{c}));
            fileNames = dir(inputDir);
            
            for f=1:length(fileNames)
                try
                    data = load(fullfile(resultDir,target,'Optimal',fileNames(f).name));
                    [val, id] = max(data.accy);
                    
                    if val > accy(o,c,chnls),
                        accy(o,c,chnls) = val;
                        
                        loc = strfind(fileNames(f).name,'_');
                        imageFile = fullfile(parentPath,'Images',target,'Optimal',[fileNames(f).name(1:loc(end)-1) '.mat']);
                        imageData = load(imageFile);
                        
                        curves{o,c,chnls} = imageData.X;
                    end
                catch
                    fprintf('Error reading file.\n');
                end
                
            end
            
        end
    end
end

%% Plot data

wave = 400:4:800;
lw = 2;

refClassifier = 'SVM';
refCounts = [1, 3, 5, 7, 9];

for o=1:length(optimalCameras)
    figure;
    subplot(2,5,1:5);
    hold on; grid on; box on;
    plot(squeeze(accy(o,:,:))','lineWidth',lw);
    xlabel('Number of channels');
    ylabel('Prediction accuracy');
    legend(classifierNames','location','southeast');
    ylim([0.2 1.05]);
    switch optimalCameras{o}
        case 'Unsup'
            title('Unsupervised');
        case 'Sup'
            title('Supervised');
    end
    
    for t=1:5
        subplot(2,5,5+t);
        hold on; grid on; box on;
        plot(wave,curves{o,cellfun(@(x) strcmp(x,refClassifier),classifierNames),refCounts(t)},'lineWidth',lw);
        xlabel('Wavelength, nm');
        title(sprintf('%i channel(s)',refCounts(t)));
        switch o
            case 1
                ylim([0 0.5]);
            case 2
                ylim([0 0.15]);
        end
    end
end

%% Plot optimal illuminants only
fs = 10;
o = 1; % Unsupervised
lw = 2;
for t=1:5
    figure;
    hold on; grid on; box on;
    plot(wave,curves{o,cellfun(@(x) strcmp(x,refClassifier),classifierNames),refCounts(t)},'lineWidth',lw);
    ylim([0 0.5]);
    set(gca,'FontSize',fs+2);
    set(gcf,'PaperPosition',[1 1 18 9]);
    set(gca,'YTick',0:0.1:0.5);
    set(gca,'XTick',400:100:800);
    
    if ~isempty(destPath)
        fName = fullfile(destPath,sprintf('Unsup_%i.eps',refCounts(t)));
        print('-depsc',fName);
    end
end


o = 2; % Supervised
lw = 2;
for t=1:5
    figure;
    hold on; grid on; box on;
    plot(wave,curves{o,cellfun(@(x) strcmp(x,refClassifier),classifierNames),refCounts(t)},'lineWidth',lw);
    ylim([0 0.15]);
    set(gca,'FontSize',fs+2);
    set(gcf,'PaperPosition',[1 1 18 9]);
    set(gca,'YTick',0:0.03:0.15);
    set(gca,'XTick',400:100:800);
    
    if ~isempty(destPath)
        fName = fullfile(destPath,sprintf('Sup_%i.eps',refCounts(t)));
        print('-depsc',fName);
    end
end

%% Plots for print and publication
lw = 2;
fs = 12;

figure;
hold on; grid on; box on;
set(gca,'FontSize',fs);
set(gcf,'PaperPosition',[1 1 9 3.5]);
dta = squeeze(accy(1:2,1,:))'*100;
hndl = plot(dta,'lineWidth',lw);
hndl(1).LineStyle = '--';
plot(1:10,ones(10,1)*85.6,'b','LineWidth',lw);
xlabel('Number of illuminants');
ylabel('Classification accuracy, %');
legend(hndl,{'Unsupervised','Supervised'},'location','southeast');
hd = annotation(gcf,'textbox',...
    [0.7 0.42 0.4 0.2],...
    'String',{'Best RGB camera'},...
    'FitBoxToText','off','lineStyle','none','fontSize',fs,'color','b');

if ~isempty(destPath)
    fName = fullfile(destPath,'NChannels.eps');
    print('-depsc',fName);
end

%% Plots for dissertation
%  This is the plot formatting I used in my dissertation.

destPath = fullfile('~','Desktop','Figures','Optimal illuminant');
if exist(destPath,'dir') == 0
    mkdir(destPath);
end

pos = [1 1 10 5];
fs = 8;
lw = 1;
ms = 5;

figure;
hold on; grid on; box on;
set(gcf,'PaperUnits','centimeters');
set(gca,'FontSize',fs);
set(gcf,'PaperPosition',pos);
set(gca','TickLabelInterpreter','LaTeX');
dta = squeeze(accy(1:2,1,:))'*100;
hndl = plot(dta,'lineWidth',lw);
hndl(1).LineStyle = ':';
hndl(1).Color = 'green';
hndl(2).Color = 'red';
plot(1:10,ones(10,1)*85.6,'b','LineWidth',lw);
xlabel('Number of illuminants','Interpreter','LaTeX','FontSize',fs);
ylabel('Classification accuracy, \%','Interpreter','LaTeX','FontSize',fs);
legend(hndl,{'Unsupervised','Supervised'},'location','southeast','Interpreter','LaTeX','FontSize',fs-2);
hd = annotation(gcf,'textbox',...
    [0.7 0.45 0.4 0.2],...
    'String',{'Best RGB camera'},...
    'FitBoxToText','off','lineStyle','none','fontSize',fs-2,'color','b',...
    'Interpreter','LaTeX');
xlim([1 size(dta,1)]);


if ~isempty(destPath)
    fName = fullfile(destPath,'NChannels.eps');
    print('-depsc',fName);
end



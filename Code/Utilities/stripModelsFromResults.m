% This is a helper script that reads in all the .mat result files, and
% removes any machine learning models that have been trained. These models
% account for the majority of the file size. The stripped files are stored
% in a separate directory.
% 
% Copyright, Henryk Blasinski 2017.

function stripModelsFromResults()

    close all;
    clear all;
    clc;

    [~, parentPath] = olRootPath();

    sourceDir = fullfile('/','Volumes','MyPassport','OptimalLight','Results');
    destDir = fullfile('/','Volumes','MyPassport','OptimalLight','ResultsSimplified');

    iterativeParse(sourceDir,destDir);
    
end


function iterativeParse(source,dest)

if exist(source,'dir') == 7
    % It is a directory, keep parsing
    contents = dir(source);
    % Remove any directory with a name starting with '.'
    contents = contents(cellfun(@(x) x(1) ~= '.',{contents(:).name}));
    
    for i=1:length(contents)
        newSource = fullfile(source,contents(i).name);
        newDest = fullfile(dest,contents(i).name);
        iterativeParse(newSource,newDest);
    end

else
    % The source is a file
    % Get it's name, load data, and save only some of it. 
    
    fprintf('Processing file: %s\n',source);
    [dirName, fName, ~] = fileparts(dest);
    loc = strfind(fName,'_');
    
    if length(loc) == 3 && isempty(strfind(source,'RedAndGreenApples'))
        if (str2double(fName(loc(2)+1:loc(3)-1)) ~= 3)
            return
        end
    end
    
    
    if exist(dirName,'dir') == 0
        mkdir(dirName);
    end
    
    load(source);
    model = cellfun(@(x) x.Trained{1},model,'UniformOutput',false);
    save(dest,'accy','predicted','true','reg','model');
end

end

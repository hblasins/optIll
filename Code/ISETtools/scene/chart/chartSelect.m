function [mRGB, mLocs, pSize, cornerPoints, mccRectHandles] = ...
    chartSelect(obj,showSelection,fullData,height,width,cornerPoints)
%
%
% Copyright Henryk Blasinski, 2014.

if ieNotDefined('obj'), obj = vcGetObject('vcimage'); end
if ieNotDefined('showSelection'), showSelection = 1; mccRectHandles = []; end
if ieNotDefined('fullData'), fullData = 0; end
if ieNotDefined('cornerPoints'), queryUser = true; else queryUser = false;end
if ieNotDefined('height'), height = 4; end
if ieNotDefined('width'), width = 6; end

% obj is either a vcimage or a sensor image
switch lower(obj.type)
    case 'vcimage'
        handles = ieSessionGet('vcimagehandles');
        dataType = 'result';
        obj = imageSet(obj,'mccRectHandles',[]);
        vcReplaceObject(obj);
        if ieNotDefined('cornerPoints')
            cornerPoints = imageGet(obj,'mcc corner points');
        end
        
    case {'isa','sensor'}
        handles = ieSessionGet('sensorWindowHandles');
        dataType = 'dvorvolts';
        obj = sensorSet(obj,'mccRectHandles',[]);
        vcReplaceObject(obj);
        if ieNotDefined('cornerPoints')
            cornerPoints = sensorGet(obj,'mcc corner points');
        end
        
    otherwise
        error('Unknown object type');
end

% If the user didn't send in any corner points, and there aren't in the
% structure, go get them from the user in the window.
if isempty(cornerPoints)
    cornerPoints = vcPointSelect(obj,4,...
        'Select (1) white, (2) black, (3) blue, (4)brown');
end

% We have cornerpoints for sure now.  Set them and draw the Rects.
switch vcEquivalentObjtype(obj.type)
    case 'VCIMAGE'
        obj = imageSet(obj,'mcc corner points',cornerPoints);
    case 'ISA'
        obj = sensorSet(obj,'mcc corner points',cornerPoints);
end
%

% Ask the user if a change is desired.  The olds one from the structure may
% not be satisfactory.
if queryUser,
    chartDrawRects(obj,'on',height,width);
    b = ieReadBoolean('Are these rects OK?');
else
    b = true;
end

if isempty(b)
    fprintf('%s: user canceled\n',mfilename);
    mRGB=[]; mLocs=[]; pSize=[]; cornerPoints=[]; mccRectHandles =[];
    return;
elseif ~b  % False, a change is desired
    switch vcEquivalentObjtype(obj.type)
        case {'VCIMAGE'}
            vcimageWindow;
        case {'ISA'};
            sensorImageWindow;
        otherwise
            error('Unknown type %s\n',obj.type);
    end
    
    % These appear to come back as (x,y),(col,row).  The upper left of the
    % image is (1,1).
    cornerPoints = vcPointSelect(obj,4,...
        'Select (1) white, (2) black, (3) blue, (4)brown');
    % should be an imageSet
    obj = imageSet(obj,'mcc corner points',cornerPoints);
end

ieInWindowMessage('',handles);

% Find rect midpoints and patch size.  mLocs are the 24 patch locations in
% (row,col) format.
[mLocs,delta,pSize] = chartRectangles(cornerPoints,height,width);

% Get the mean RGB data or the full data from the patches in a cell array
% The processor window is assumed to store linear RGB values, not gamma
% corrected.
mRGB = chartPatchData(obj,mLocs,delta,fullData,dataType,height,width);

% Plot the rectangle that encloses these points.
if showSelection, chartDrawRects(obj,[],height,width); end

ieInWindowMessage('',handles);

return;


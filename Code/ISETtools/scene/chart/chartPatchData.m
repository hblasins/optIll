function mRGB = chartPatchData(obj,mLocs,delta,fullData,dataType,height,width)
%
% Copyright, Henryk Blasinski, 2014
%

if ieNotDefined('obj'),   error('vcimage or sensor required'); end
if ieNotDefined('mLocs'), error('Mid locations required'); end
if ieNotDefined('delta'), error('Patch spacing required'); end
if ieNotDefined('fullData'),fullData = 0; end         % Mean, not all the points
if ieNotDefined('delta'),   dataType = 'result'; end  % Default for vcimage
if ieNotDefined('height'), height = 4; end
if ieNotDefined('width'), width = 6; end

if fullData  % Every value in the patch
    mRGB = cell(1,height*width);
    for ii = 1:height*width
        % mLocs(:,mPatch) is a column vector with (row,col)' for the
        % mPatch.
        theseLocs = macbethROIs(mLocs(:,ii),delta);
        mRGB{ii} = vcGetROIData(obj,theseLocs,dataType);
    end
else  % Mean values from each patch
    mRGB = zeros(height*width,3);
    for ii = 1:height*width
        % mLocs(:,mPatch) is a column vector with (row,col)' for the
        % mPatch.
        % This code doesn't work properly for the case of an image sensor.
        % It needs to protect against the NaNs returned in that case.  It
        % works OK for the vcimage.  Fix this some day.
        theseLocs = macbethROIs(mLocs(:,ii),delta);
        mRGB(ii,:) = mean(vcGetROIData(obj,theseLocs,dataType));
    end
end

return

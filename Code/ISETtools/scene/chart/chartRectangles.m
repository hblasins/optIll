function [mLocs,delta,pSize] = chartRectangles(cornerPoints,height,width)
%
% Copyright Henryk Blasinski 2014.

if ieNotDefined('cornerPoints'), error('Point corners required'); end
if ieNotDefined('height'), height = 4; end
if ieNotDefined('width'), width = 6; end

% cornerpoints are (x,y) (col,row) format.
cornerPoints = fliplr(cornerPoints);
 
% cornerPoints contains the positions of the location of the four
% corners.
mWhite = cornerPoints(1,:);
mBlack = cornerPoints(2,:);
mBlue  = cornerPoints(3,:);
mBrown = cornerPoints(4,:);

% Find the affine transformation that maps the selected point values into a
% standard spatial coordinate frame for the Macbeth Color Checker, with the
% white in the lower left, black in lower right, and brown on the upper
% left.
offset = mWhite;
mBlack = mBlack - offset;
mBlue  = mBlue  - offset;
mBrown = mBrown - offset;

% Find the linear transformation that maps the non-white points into the
% ideal positions.  These are (x,y) format.
%  White -> 0,0
%  Black -> width,0
%  Blue ->  width,height
%  Brown -> 0,height
ideal = [width width 0; 0 height height];
current = [mBlack(:), mBlue(:), mBrown(:)];

%  current = L * ideal
L = current*pinv(ideal);

%  Any coordinate in the ideal target can be transformed into the current
%  row, col values in the current data by currentLoc = L*idealLoc
%  So, for example, the red is at 2.5,1.5
%  In the current data this would be
% (L*[2.5,1.5]') + offset(:)
%
% So now, we make up the coordinates of all 24 patches.  These are
[X,Y] = meshgrid((0.5:1:(width-0.5)),(0.5:1:(height - 0.5)));
idealLocs = [X(:),Y(:)]';

% The mLocs contains (rows,cols) of the center of the 24 patches.
% These are from white (lower left) reading up the first row, and then back
% down to the 2nd column, starting at the gray, and reading up again
mLocs = round(L*idealLocs + repmat(offset(:),1,height*width));
flipIt = 1:height*width;
flipIt = reshape(flipIt,height,width);
flipIt = flipud(flipIt);
flipIt = flipIt(:);
mLocs = mLocs(:,flipIt);

% Build a square of a certain size around the point mLocs(:,N)
% We need to know whether the white->black is down the rows (first
% dimension) or down the columns (second dimension). If the row difference
% between the first two is much larger than the column difference, then we
% assume white black is down the rows.  Otherwise, we assume white->black
% is across the columns.

deltaX = round(abs(cornerPoints(1,2) - cornerPoints(2,2))/width);
deltaY = round(abs(cornerPoints(1,1) - cornerPoints(4,2))/height);


% We want to pick out a square region that is within the size of the patch.
% So, we divide the estimated width by 3 and take the smaller one.  This
% way, the coverage is about 2/3rds of the estimated square size.  If we
% divide by 4, rather than 3, we will get 1/2 the patch size.
delta = round(min(deltaX,deltaY)/3);
pSize = 2*delta + 1;

% Debug:
% Put up the mean locations in the sensor image
% plot(mLocs(2,:),mLocs(1,:),'wo')
  
return;

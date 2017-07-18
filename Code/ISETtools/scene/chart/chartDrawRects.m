function obj = chartDrawRects(obj,onoff,height,width)


if ieNotDefined('obj'), error('Structure required'); end
if ieNotDefined('onoff'), onoff = 'on'; end  % Default is on
if ieNotDefined('height'), height = 4; end
if ieNotDefined('width'), width = 6; end

switch onoff
    case 'on'
        % Should be an imageGet(obj,'mcc corners');
        switch vcEquivalentObjtype(obj.type)
            case 'VCIMAGE'
                cornerPoints = imageGet(obj,'mcc corner points');
                a = get(vcimageWindow,'CurrentAxes');
                
            case 'ISA'
                cornerPoints = sensorGet(obj,'mcc corner points');
                a = get(sensorImageWindow,'CurrentAxes');
                
            otherwise
                error('Unknown object type %s',obj.type);
        end
        if isempty(cornerPoints), error('No mcc corner points'); end
        
        [mLocs, delta] = chartRectangles(cornerPoints,height,width);
        
        % Plot the rectangles
        rectHandles = zeros(height*width,1);
        for ii=1:height*width
            theseLocs = macbethROIs(mLocs(:,ii),delta);
            corners = convhull(theseLocs(:,1),theseLocs(:,2));
            hold(a,'on');
            rectHandles(ii) = plot(a,theseLocs(corners,2),theseLocs(corners,1),...
                'Color',[1 1 1], 'LineWidth',2);
        end
        
        % Store rectHandles so we can delete them later.
        % A refresh will delete them also, apparently.
        switch lower(obj.type)
            case 'vcimage'
                obj = imageSet(obj,'mccRectHandles',rectHandles);
            case {'isa','sensor'}
                obj = sensorSet(obj,'mccRectHandles',rectHandles);
        end
        vcReplaceObject(obj);
    
    case 'off'
        % Delete handles from current axis and update the object
        switch lower(obj.type)
            case 'vcimage'
                rects = imageGet(obj,'mcc Rect Handles');          
                delete(rects);
                obj = imageSet(obj,'mccRectHandles',[]);

            case {'isa','sensor'}
                rects = sensorGet(obj,'mcc Rect Handles');
                delete(rects);
                obj = sensorSet(obj,'mccRectHandles',[]);
        end
        vcReplaceObject(obj);

    otherwise
        error('Unknown on/off %s\n',onoff);
end

return;

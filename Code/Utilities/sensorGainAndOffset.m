function [gain, offset] = sensorGainAndOffset(scene,oi,sensor)

% [ gain, offset ] = sensorGainAndOffset(scene, oi, sensor)
%
% This function computes the gain and offset that relate the linear image
% formation model with pixel intensities. In general a pixel intensity can
% be expressed as a linear model
%
% p = gain*[(reflectance*illuminant*responsivity)*deltaWavelength] + offset
%
% where p is the pixel intensity. The quantity in the square brackets is
% the sum over wavebands of the product of spectral reflectance, the
% illuminant spectral power distribution and the sensor responsivity.
% The offset can be due to DSNU, noise or analog offset introduced by the
% sensor.
%
% We assume that p is normalized in [0, 1], so that the computed gain
% and offset values are equally valid for normalized voltages or digital 
% values.
% 
% Copyright, VISTA Lab 2016

% Compute the conversion from scene radiance to irradaince, this depends on
% the scene distance, fNumber and magnification.
sDist = sceneGet(scene,'distance');
fN    = opticsGet(oiGet(oi,'optics'),'fNumber');     
m     = opticsGet(oiGet(oi,'optics'),'magnification',sDist);
rad2irr = pi /(1 + 4*fN^2*(1+abs(m))^2);




% Estimate how many points in the scene contributes to a single pixel.
gridSpacing = 1/sensorGet(sensor,'nSamplesPerPixel');
gridSpacing = 1/round(1/gridSpacing);
nGridSamples = 1/gridSpacing;

if nGridSamples == 1, pdArray = pixelGet(sensorGet(sensor,'pixel'),'fillfactor');
else                  pdArray = sensorPDArray(sensor,gridSpacing);
end
area = pixelGet(sensorGet(sensor,'pixel'),'area')/(nGridSamples^2);

% Convert photo-electrons to voltage.
pixel = sensorGet(sensor,'pixel');
q = vcConstants('q');

irr2curr = 1/q;
cur2volt = pixelGet(pixel,'conversionGain') * sensorGet(sensor,'integrationTime');
volt2dv = 1/pixelGet(pixel,'voltageswing');
analogGain = sensorGet(sensor,'analoggain');


gain = q*volt2dv*cur2volt*irr2curr*rad2irr*pdArray*area/analogGain;


%% Offset computations depend on the noise model used
noiseFlag = sensorGet(sensor,'noiseFlag');

offset = 0;
if noiseFlag >=0
    % Account for the analog offset only
    offset = offset + sensorGet(sensor,'analogoffset')/sensorGet(sensor,'analoggain');
end
if noiseFlag >=2
    % We need to include the effects of dark voltage
    offset = offset + pixelGet(pixel,'dark Voltage')*sensorGet(sensor,'integrationTime')/sensorGet(sensor,'analoggain');
end

% Convert the offset to the [0 1] scale.
offset = offset*volt2dv;
end
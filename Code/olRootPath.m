function [rootPath, parentDir] = olRootPath()

% function [rootPath, parentDir] = olRootPath()
%
% Returns the absolute path for the directory containing the optimal
% illuminant computation code.
%
% Copyright, Henryk Blasinski 2017

rootPath = which('olRootPath');
rootPath = fileparts(rootPath);

id = strfind(rootPath,'/');
parentDir = rootPath(1:(id(end)-1));

return

function [scStruct] = scoreSCdataset(scStruct, metModel)
% Compute RAS with a scStruct as input
% 
% USAGE:
%
%   [scStruct] = scoreSCdataset(scStruct, metModel)
%
%
% INPUT:
%   scStruct:   Single Cells dataset in a structure built with makeSCdataset function.
%   metModel:   Name of field corresponding to the model in scStruct upon which integrate RAS
%
% OUTPUTS:
%   scStruct:   Single Cells dataset in a structure with RAS field
%
%
% .. Author:
%       - Davide Maspero 30/01/2018

try % parallel mode
    scStruct.RAS = getDataScored(metModel, scStruct.GenName, scStruct.TPMsc', scStruct.CellType,  1);
catch
    scStruct.RAS = getDataScored(metModel, scStruct.GenName, scStruct.TPMsc', scStruct.CellType,  0);
end

end
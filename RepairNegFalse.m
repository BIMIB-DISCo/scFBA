function [scStruct] = RepairNegFalse(scStruct, epsilon)
% Preprocess of the expression level to handle the noise in the single
% cells data.
% 
% USAGE:
%
%   [scStruct] = RepairNegFalse(scStruct, epsilon)
%
% INPUT:
%   scStruct:       Single Cells dataset in a structure built with makeSCdataset function.
%   epsilon:        Value added to single cell transcript levels.
%
% OUTPUTS:
%   scStruct:       Single Cells dataset in a structure built with makeSCdataset function.
%
%
% .. Author:
%       - Davide Maspero 30/01/2018

if nargin < 2
    epsilon = min(min(scStruct.TPMsc));
end

if ~isfield(scStruct, 'TPMpl')
    error('No Pooled or Bulk data founded. scStruct must have TPMpl field.')
end

count = 0;
%level 1
for i=1:length(scStruct.TPMpl)
    if((scStruct.TPMpl(i) > epsilon) && (mean(scStruct.TPMsc(:,i)) < epsilon*2))
        scStruct.TPMsc(:,i) = scStruct.TPMpl(i);
        count = count + 1;
    end
end
disp(strcat(num2str(count),{' gene(s)'}, {' modified'}));

end
function [scStruct] = makeSCdataset(PoolExp, ExpGxC_SC, CellType, Gene_Names, epsilon)
% Preprocess dataset of single cell transcriptomic profiles to run scFBA
% functions.
% 
% USAGE:
%
%   scStruct = makeSCdataset(PoolExp, ExpGxC_SC, CellType, Gene_Names, epsilon)
%
%
% INPUT:
%   PoolExp:            Vector with expression profile of Pooled or Bulk
%                       cells.
%   ExpGxC_SC:          matrix genes x cells with TPM or RPKM or any other
%                       values for the single cell expression levels NOT
%                       log scaled.
%   CellType:           Identifier of each single cells.
%   Gene_Names:         Identifier of each genes. The format must be the
%                       same as the gene identifier used in metabolic model.
% OPTIONAL INPUTS:
%   epsilon:            value to add on each single cells expression level
%
% OUTPUTS:
%   scStruct:           Single Cells dataset in a structure built with makeSCdataset function.


if nargin < 5
    epsilon = 0;
end

if size(CellType, 1)>1
    scStruct.CellType = CellType;
else
    scStruct.CellType = CellType';
end
if size(Gene_Names, 1)>1
    scStruct.GenName = Gene_Names;
else
    scStruct.GenName = Gene_Names';
end
if size(PoolExp, 1)==1
    scStruct.TPMpl = PoolExp;
else
    scStruct.TPMpl = PoolExp';
end
if size(ExpGxC_SC, 1)==size(CellType, 1)
    scStruct.TPMsc = ExpGxC_SC;
else
    scStruct.TPMsc = ExpGxC_SC';
end

scStruct.TPMsc = scStruct.TPMsc + epsilon;

end



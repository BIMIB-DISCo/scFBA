function [datasetOut, datasetCell, listGenNotInModel, duplicateGenes] = extractGeneRow(model, datasetIn, datasetCol, dict)
% Filter out gene not present in metabolic model from a dataset,
%
% USAGE:
%
%       [datasetOut, datasetCell] = extractGeneRow(model, datasetIn, datasetCol, dict)
%
% INPUT:
%   model:          metabolic model in COBRA format
%   datasetIn:      matrix genes x samples in any MATLAB table or cell arrays format
%   datasetCol:     number of column with gene IDs to match gene IDs in
%                   metabolic model or in dictionary    
%
%
% OPTIONAL INPUTS:
%   dict:           cell arrays with first column contain dataset gene IDs,
%                   second column contain model gene IDs
%
%
% OUTPUTS:
%   datasetOut:     Dataset filtered on metabolic genes.
%   datasetCell:    Dataset filtered in cell array format.
%
%
% .. Author:
%       - Davide Maspero 30/01/2018


conv = false;
if istable(datasetIn)
    if isnumeric(datasetIn.(datasetCol))
        datasetIn.(datasetCol) = strtrim(cellstr(num2str(datasetIn.(datasetCol))));
    end
    varName = datasetIn.Properties.VariableNames;
    datasetIn = table2cell(datasetIn);
    datasetIn(:,datasetCol) = cellstr(datasetIn(:,datasetCol));
    conv = true;
end
if nargin > 3 %if dict was pass
    if istable(dict)
        if isnumeric(dict.(1))
            dict.(1) = strtrim(cellstr(num2str(dict.(1))));
        end
        idConvName = dict.Properties.VariableNames;
        dict = table2cell(dict);
    else
        idConvName = [{'OrigID'}, {'ConvID'}];
    end
    dict = cellstr(dict);
    dict = dict(ismember(dict(:,2), model.genes),:);
    if size(dict, 1) ~= size(model.genes, 1)
        error('Gene ID dictionary is not compleate, please check');
    end
else
    dict = model.genes;
end

% [ismem, idxDataDict] = ismember(dict(idxDictModel,1), datasetIn(:,datasetCol));
% idxDataDict=idxDataDict(ismem);
% idxDictModel=idxDictModel(ismem);
% datasetOut = datasetIn(idxDataDict,:);
[ismem, idxDataDict] = ismember(datasetIn(:,datasetCol), dict(:,1));
idxDictModel=idxDataDict(idxDataDict~=0);
datasetOut = datasetIn(ismem,:);
%datasetOut = datasetIn(ismember(datasetIn(:,datasetCol), dict(idxDictModel,1)),:);

if nargin > 3
    genNotinModel = ismember(dict(:,1), datasetOut(:,datasetCol))==0;
    genNotinModel = find(ismember(model.genes, dict(genNotinModel,2)));
else
    genNotinModel = find(ismember(model.genes, datasetOut(:,datasetCol))==0);
end
%genNotinModel = find(genNotinModel==0);

if ~isempty(genNotinModel)
    disp(strcat(num2str(length(genNotinModel)), ' gene(s) in model but not in dataset:'))
    listGenNotInModel = model.genes(genNotinModel);
    disp(listGenNotInModel)
end

duplicateGenes = datasetOut(:,datasetCol);
uni = unique(duplicateGenes);
[~, idxUni] = ismember(uni, duplicateGenes);
duplicateGenes(idxUni) = [];

if ~isempty(duplicateGenes)
    disp(strcat(num2str(length(duplicateGenes)), ' gene(s) duplicate in dataset:'));
    disp(duplicateGenes);
end

if conv
    datasetOut = cell2table(datasetOut,'VariableNames',varName);
    if nargin > 3
        dict = cell2table(dict(:,2),'VariableNames',idConvName(2));
        datasetOut = [dict(idxDictModel,1) datasetOut];
    end
else
    if nargin > 3
        datasetOut = [dict(idxDictModel) datasetOut];
    end
end

warning('off', 'MATLAB:table:RowsAddedExistingVars');

if ~isempty(genNotinModel)
    datasetOut(size(datasetOut, 1)+1 : size(datasetOut, 1)+ size(genNotinModel, 1),1) = model.genes(genNotinModel);
    for i=1:size(datasetOut, 2)
        if isnumeric(datasetOut{end,i})
            IdxNan(i) = 1;
        end
    end
    datasetOut{end-size(genNotinModel, 1)+1:end,(IdxNan==1)} = NaN;
    IdxNan(1)=1;
    datasetOut{end-size(genNotinModel, 1)+1:end,(IdxNan==0)} = {'missing'};
end
if nargout > 1
    ColName = datasetOut.Properties.VariableNames;
    datasetCell = table2cell(datasetOut);
    datasetCell = [ColName; datasetCell];
end
warning('on', 'MATLAB:table:RowsAddedExistingVars');
end
function [scStructOut] = single2IntPopModel(scStruct, singleModel, IdxExRxns, IdxCoopRxn, compExt)
% Build a scStruct with modelPop, modelDel, modelFVA and modelFVAInt from a
% metabolic model in cobra format.
% 
% USAGE:
%
%   [scStructOut] = single2IntPopModel(scStruct, singleMetModel, IdxExRxns, IdxCoopRxn, compExt)
%
%
% INPUT:
%   scStruct:           Single Cells dataset in a structure built with makeSCdataset function.
%   singleMetModel:     Metabolic model in COBRA format
%   IdxExRxns:          ID reactions to be cosidered es exchange in
%                       modelPop. For details see createPopModel help or
%                       popFBA paper (ha senso?)
%   IdxCoopRxn:         ID reactions to be cosidered es cooperative in
%                       modelPop. For details see createPopModel help or
%                       popFBA paper (ha senso?)
%   compExt:            letter of the shared compartment of the population in modelPop. 
%                       For details see createPopModel help or popFBA paper (ha senso?)
%
% OUTPUTS:
%   scStructOut:        struct with following fields added.
%                           'modelPop':     Single model multiplied as many 
%                                           times as the cells in the scStruct 
%                                           (lenght of cellType field)
%                           'modelDel':     (optional) modelPop with scStruct.GenNameOff
%                                           deleted, if that field exist.
%                                           GenNameOff is build by Genes_Sign function.
%                           'modelFVA':     ModelDel or ModelPop with boundaries setted up
%                                           by fluxVariability
%                           'modelFVAInt':  modelFVA with boundaries
%                                           multiplied by RAS and
%                                           normalized over the sum of RASs
%                                           of the whole population.
%                           'RAS':          table reactions x cells with
%                                           RAS computed

numSC = length(scStruct.CellType);
% exchange boundaries multiply for the numbers of single cells
singleModel.lb(IdxExRxns) = singleModel.lb(IdxExRxns)*numSC;
singleModel.ub(IdxExRxns) = singleModel.ub(IdxExRxns)*numSC;

% delete gene indicate in struct.genNameOff
if ~isfield(singleModel, 'rxnGeneMat')
    try
        singleModel = buildRxnGeneMat(singleModel);
    catch
        error('Matrix with reactions - genes assotiation (rxnGeneMat) not found in singleModel and it is impossible to build. Please check.')
    end
end

% building popModel
disp('Building population model');
scStruct.modelPop = createPopModel(singleModel, IdxExRxns, IdxCoopRxn, numSC, compExt);

deleted = false;

if isfield(scStruct, 'GenNameOff')
    if ~isempty(scStruct.GenNameOff)
        deleted = true;
        disp('Delete gene(s) in population model');
        scStruct.modelDel = deleteModelGenes(scStruct.modelPop, scStruct.GenNameOff);
    end
end

% compute RAS over the modelDel

if isfield(scStruct, 'RAS')
    if isempty(scStruct.RAS)
        disp('Compute RAS')
        scStruct = scoreSCdataset(scStruct, singleModel);
    end
else
    disp('Compute RAS')
    scStruct = scoreSCdataset(scStruct, singleModel);
end

disp('Integrate RAS in model');
if deleted
    scStruct = integrateRAS(scStruct, 'modelDel');
else
    scStruct = integrateRAS(scStruct, 'modelPop');
end

scStructOut = scStruct;

end


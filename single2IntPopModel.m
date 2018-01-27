function [scStructOut] = single2IntPopModel(scStruct, singleMetModel, IdxExRxns, IdxCoopRxn, compExt)

numSC = length(scStruct.CellType);
% exchange boundaries multiply for the numbers of single cells
singleMetModel.lb(IdxExRxns) = singleMetModel.lb(IdxExRxns)*numSC;
singleMetModel.ub(IdxExRxns) = singleMetModel.ub(IdxExRxns)*numSC;

% building popModel
disp('Building population model');
scStruct.modelPop = createPopModel(singleMetModel, IdxExRxns, IdxCoopRxn, numSC, compExt);

% delete gene indicate in struct.genNameOff
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
        scStruct = scoreSCdataset(scStruct, singleMetModel);
    end
else
    disp('Compute RAS')
    scStruct = scoreSCdataset(scStruct, singleMetModel);
end

disp('Integrate RAS in model');
if deleted
    scStruct = integrateRAS(scStruct, 'modelDel');
else
    scStruct = integrateRAS(scStruct, 'modelPop');
end

scStructOut = scStruct;

end


function [scStructOUT, vectRAS] = integrateRAS(scStruct, fieldName)
% Constrain model by multiplying each reaction boundaries times theirs RAS.
% RAS must be computed with scoreSCdataset function.
%
% USAGE:
%
%   [scStructOUT, vectScore] = integrateRAS(scStruct, fieldName)
%
% INPUT:
%   scStruct:       Single Cells dataset in a structure built with makeSCdataset function.
%   fieldName:      Name of field corresponding to the model in scStruct upon which integrate RAS
%
%
% OUTPUTS:
%   scStructOUT:    cell arrays or MATLAB table converted.
%   vectRAS:        vector with RAS compute for each reactions in model

modelPop = getfield(scStruct, fieldName);
nPop = length(scStruct.CellType);

IdxPop_0 = find(endsWith(modelPop.rxns, '_0')==1);
RxnList = modelPop.rxns(IdxPop_0);

%compute fluxVariability and save Fmin, FMax vectors
objVect = modelPop.c;
modelPop.c(:) = 0;
[Fmin,FMax] = fluxVariability(modelPop, 100, '', RxnList); %controllare che 100 vada bene
modelPop.c = objVect;

Bound = max(abs(Fmin), abs(FMax));

modelPopBound = modelPop;
%
%remove _nPop
for k=1:length(RxnList)
    while (RxnList{k}(end) ~= '_')
        RxnList{k} = RxnList{k}(1:end-1);
    end
end

for j=0:nPop-1
    tmpRxnList = RxnList;
    tmpRxnList = strcat(tmpRxnList, num2str(j));
    for i=1:length(tmpRxnList)
        if modelPop.lb(IdxPop_0(i))<0
            if modelPop.ub(IdxPop_0(i)) > 0 %rev
                modelPopBound = EditBoundaries(modelPopBound, tmpRxnList(i), -Bound(i), Bound(i), false, true);
            else % Solo Bw
                modelPopBound = EditBoundaries(modelPopBound, tmpRxnList(i),  -abs(Fmin(i)), 0, false, true); % mi assicuro che il lb sia negativo
            end
        else %Solo Fw
            modelPopBound = EditBoundaries(modelPopBound, tmpRxnList(i), 0, abs(FMax(i)), false, true); % mi assicuro che l'ub sia positivo
        end
    end
end

if isfield(scStruct, 'RAS')
    if isempty(scStruct.RAS)
        error('RAS not found in single cell struct. Please use the function scoreSCdataset before run this function.');
        return
    end
else
    error('RAS not found in single cell struct. Please use the function scoreSCdataset before run this function.');
    return
end
scStructOUT = scStruct;
scStructOUT.modelFVA = modelPopBound;

modelInt = modelPopBound;
modelBck = modelPopBound;
FM = max(modelInt.ub);
Fm = min(modelInt.lb);

vectRAS = zeros(length(modelBck.rxns), 1);

for i=1:length(modelBck.rxns)
    tmpSpl = split(modelBck.rxns(i), '_');
    if(~strcmp(tmpSpl(end), '#'))
        if(strcmp(modelBck.rules(i),''))
            vectRAS(i) = 0;
        else
            %rules = modelBck.grRules{i};
            %MaxCoef = getObjective(rules, dataset.GenName, maxTransc);
            Coeff = scStruct.RAS{modelBck.RxnsOldIdx(i), str2double(tmpSpl{end})+2};
            %somma tutti i coefficenti dati dalla risoluzione della regola per ogni cellula
            sumCoeff = sum(scStruct.RAS{modelBck.RxnsOldIdx(i), 2:end});
            vectRAS(i) = Coeff / sumCoeff;
            
            if isnan(vectRAS(i))
                continue
            end
            
            if(modelInt.lb(i) > Fm && abs(modelInt.lb(i)) > 10^-3) %se è = a infinito o inferiore al minimo non modifico i vincoli.
                %modelInt.lb(i) = modelInt.lb(i) * vectRAS(i);
                modelInt.lb(i) = mapInRange(vectRAS(i), 0,1, -10^-3, modelInt.lb(i));
            end
            if(modelInt.ub(i) < FM && abs(modelInt.ub(i)) > 10^-3)
                %modelInt.ub(i) = modelInt.ub(i) * vectRAS(i);
                modelInt.ub(i) = mapInRange(vectRAS(i), 0,1, 10^-3, modelInt.ub(i));
            end
        end
    end
end

scStructOUT.modelFVAInt = modelInt;

end
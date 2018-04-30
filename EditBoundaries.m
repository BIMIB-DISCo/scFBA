function [modelOut, TableRes] = EditBoundaries(model,rxn, newLb, newUb, NotSure, exactMatch)
% Divide the popFBA flux solution in a matrix nReaction x nPop 
%
% USAGE:
%
%   [modelOut, TableRes] = EditBoundaries(model,rxn, newLb, newUb, NotSure, exactMatch)
%
% INPUT:
%   model:          metabolic model in COBRA format
%   rxn:            sub string to search in model.rxns field
%
% OPTIONAL INPUTS:
%   newLb:          new lower bound for the reactions founded, if empy values ([]) is
%                   passed no change occours for this bound
%   newUb:          new upper bound for the reactions founded, if empy values ([]) is
%                   passed no change occours for this bound
%   NotSure:        if TRUE a message box appeare before change anything (Default =  TRUE)
%   exactMatch:     if TRUE the string rxn is compeared with the complete model.rxns strings
%                   with FALSE the match is made between rxn and any model.rxns substrings 
%                   (Default = FALSE)
%
%
% OUTPUTS:
%   modelOut:       New model in COBRA format with boundaries changed.
%   TableRes:       A table with id, reaction identifier and boundaries
%                   founded. Usefull to keep the vector position of the
%                   model.rxns founded.
%
%
% .. Author:
%       - Davide Maspero 30/01/2018


if nargin < 6
    exactMatch = false;
end
if exactMatch
    idx = find(strcmp(model.rxns, rxn)==1);
else
    idxS =  strfind(lower(model.rxns),lower(rxn));
    idx = find(not(cellfun('isempty', idxS)));
end
%TableRes = table();
if nargin < 3
    if(isempty(idx))
        disp(rxn);
        disp('No reaction found');
        modelOut = model;
        return
    else
        TableRes = table(idx, model.rxns(idx), model.lb(idx),  model.ub(idx), 'VariableNames', {'ID', 'Reaction', 'LowerBound', 'UpperBound'})
        disp('No change was made');
        modelOut = model;
        return
    end
elseif nargin < 5
    NotSure = true;
end

if(isempty(idx))
    disp(rxn);
    disp('No reaction found');
    modelOut = model;
elseif(NotSure)
    TableRes = table(idx, model.rxns(idx), model.lb(idx),  model.ub(idx), 'VariableNames', {'ID', 'Reaction', 'LowerBound', 'UpperBound'})
    choice = questdlg('Would you like to edit these boundaries?', 'Edit Boundaries', 'Yes','No','No');
    switch choice
        case 'Yes'
            if ~isempty(newLb)
                model.lb(idx) = newLb;
            end
            if ~isempty(newUb)
                model.ub(idx) = newUb;
            end
            VetRev = zeros(length(idx), 1); %change the reversibility of the reactions
            VetRev(model.lb(idx) < 0 & model.ub(idx) > 0) = 1;
            model.rev(idx) = VetRev;
            modelOut = model;
        case 'No'
            modelOut = model;
    end
else
    if nargout > 1
    TableRes = table(idx, model.rxns(idx), model.lb(idx),  model.ub(idx), 'VariableNames', {'ID', 'Reaction', 'LowerBound', 'UpperBound'});
    end
    if ~isempty(newLb)
        model.lb(idx) = newLb;
    end
    if ~isempty(newUb)
        model.ub(idx) = newUb;
    end
    VetRev = zeros(length(idx), 1); %change the reversibility of the reactions
    VetRev(model.lb(idx) < 0 & model.ub(idx) > 0) = 1;
    model.rev(idx) = VetRev;
    modelOut = model;
end
end

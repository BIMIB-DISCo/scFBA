function [TableRes] = getFluxes(model, flux, rxn, hBins, exactMatch)
%User friendly function to visualize the fluxes of a FBA model or a
%population model

%flux could be either the struct result from optimizeCbModel or only a
%vector of fluxes

if ~isstruct(model)
    reactions = model;
else
    reactions = model.rxns;
end

if nargin < 4
    hBins = 0;
end
if nargin < 5
    exactMatch = false;
end
if strcmp(rxn, 'all reactions')
    idx = [1:length(reactions)]';
else
    if exactMatch
        idx = find(strcmp(reactions, rxn)==1);
    else
        idxS =  strfind(lower(reactions),lower(rxn));
        idx = find(not(cellfun('isempty', idxS)));
    end
end
TableRes = table();

if(isempty(idx))
    disp(rxn);
    disp('No reaction found');
    return
else
    if isstruct(flux)
        TableRes = table(idx, reactions(idx), flux.x(idx), 'VariableNames', {'ID' 'Reaction', 'Flux',});
    else
        TableRes = table(idx, reactions(idx), flux(idx,:), 'VariableNames', {'ID' 'Reaction', 'Flux',});
    end
    if hBins > 0
        hist(TableRes.Flux, hBins);
    end
    return
end
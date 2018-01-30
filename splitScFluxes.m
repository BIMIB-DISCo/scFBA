function [FluxPop, RxnPop] = splitScFluxes(model, flux, nPop)
% Divide a population model optimization result in a matrix reactions x
% nPop
%
% USAGE:
%
%   [scStruct] = Genes_Sign(scStruct, CutOff, Hist, ExpXLS, Epsilon, nBins, FitFunction)
%
% INPUT:
%   model:      Population model with model.rxns field.
%   flux:       Optimization of model compute with optimizeCbModel function
%
% OUTPUT:
%   FluxPop:    Matrix reactions x nPop fluxes
%   RxnPop:     Reactions identifier in the same order as the FluxPop rows.
%
%
% .. Author:
%       - Davide Maspero 30/01/2018

VRxnName = model.rxns;
if isstruct(flux)
    VFluxVal = flux.x;
else
    VFluxVal = flux;
end

%creazione matrice con indici della popolazione

FluxPop = [];

[VRxnName, IdxSort] = sort(VRxnName);
VFluxVal = VFluxVal(IdxSort);

for i=0:nPop-1
    Suffix = strcat('_', num2str(i));
    IdxPop = find(endsWith(VRxnName, Suffix)==1);
    FluxPop = [FluxPop VFluxVal(IdxPop, :)];
end

RxnPop = VRxnName(IdxPop, :);

for k=1:length(RxnPop)
    while (RxnPop{k}(end) ~= '_')
        RxnPop{k} = RxnPop{k}(1:end-1);
    end
    RxnPop{k} = RxnPop{k}(1:end-1);
end
%    cg = clustergram(FluxPop'); % cluster hitmap cellule x flussi
%     Y = pdist(FluxPop');
%     Z = linkage(Y);
%     dendrogram(Z, 0);
end












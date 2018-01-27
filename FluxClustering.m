function [FluxPop, RxnPop] = FluxClustering(model, flux, nPop)

VRxnName = [];
VFluxVal = [];
%Rimerge reazioni splittate in Fw e BW
for i=1:length(model.rxns)
    if(strcmp(model.rxns{i}(1:3), 'Fw_'))
        RxName = model.rxns{i}(4:end);
        FluxVal = flux(i,:) - flux(i+1,:); %Fw - Bw
        VRxnName = [VRxnName; RxName];
        VFluxVal = [VFluxVal; FluxVal];
    elseif(~strcmp(model.rxns{i}(1:3), 'Bw_'))
        RxName = model.rxns(i);
        FluxVal = flux(i,:);
        VRxnName = [VRxnName; RxName];
        VFluxVal = [VFluxVal; FluxVal];
    end
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












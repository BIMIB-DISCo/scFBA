%% test scFBA functions

% load transcript matrix T
load('data/BC04dataset.mat')

% load template metabolic network A
HMRcore = readCbModel('HMRcore.xml');

% Load dictionary gene id Ensemble to Hugo Name
load('data/DictCORE_ENS2HGNC.mat');

% define exchange reactions IDs
[~, Ex_id] = EditBoundaries(HMRcore, 'Ex_');
IdxExRxns = Ex_id.ID;

% define cooperation reactions IDs
[~, Coop_id] = EditBoundaries(HMRcore, '_COOP');
[~, Biomass_id] = EditBoundaries(HMRcore, 'biomass_synthesis');
IdxCoopRxn = [Coop_id.ID; Biomass_id.ID];

% filter out genes not in HMRcore model
BC04_filt = extractGeneRow(HMRcore, BC04dataset, 1, DictCORE);

% Organize single cell samples in a data structure
BC04 = makeSCdataset(BC04_filt.BC04_Pooled, BC04_filt{:,7:end}, BC04_filt.Properties.VariableNames(7:end), BC04_filt.HGNC_ID, 10^-4);

% Determine genes that have 0 transcript in both pooled and single cells
% samples (GenesOff). That genes will be complete delete from the model.
BC04 = Genes_Sign(BC04);

% Reduce the noise of single cells data by compare their expression level
% with pooled.
BC04 = RepairNegFalse(BC04);

% Initialize COBRA toolbox
% try
%     changeCobraSolver('gurobi');
% catch
%     changeCobraSolver('glpk');
% end

changeCobraSolver('glpk');

% Modify constraints as in Damiani et al
HMRcore = ScFBAExpSetting(HMRcore, length(BC04.CellType));

% Compute RAS and integrate them in modelFVAInt. During the process also
% modelPop (template model multiplied many times as single cell samples),
% modelDel (with GenesOff deleted), modelFVA (with boundaries determinated
% from fluxVariability analisys) are created and added to the scStructure
BC04 = single2IntPopModel(BC04, HMRcore, IdxExRxns, IdxCoopRxn, 's');

% integrated Model can now be optimize
fluxIntBC04 = optimizeCbModel(BC04.modelFVAInt);

% Split flux vector of the optimal solution in a matrix reaction x nPop fluxes, one column for each single
% cells.
[FluxPopBC04, RxnPop] = splitScFluxes(BC04.modelFVAInt, fluxIntBC04, length(BC04.CellType));

% normalize fluxes between 0 and 1.
FluxPopNormBC04 = FluxPopBC04;
for i=1:length(RxnPop)
    FluxPopNormBC04(i,:) = ((FluxPopBC04(i,:) - min(FluxPopBC04(i,:)))./(max(FluxPopBC04(i,:)) - min(FluxPopBC04(i,:))));
end

% Draw a clustergram of the fluxes compute
clst = clustergram(FluxPopNormBC04( ~isnan(FluxPopNormBC04(:,1)) ,:),'Cluster',2,'Colormap',redbluecmap,'symmetric',false);
addTitle(clst, 'BC04 fluxes');





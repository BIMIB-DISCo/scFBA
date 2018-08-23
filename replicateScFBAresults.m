% Script replicate scFBA paper
%% Load gene ID dictionary
load('data/DictCORE_ENS2HGNC.mat');

% Import sbml model HMRcore
HMRcore = readCbModel('HMRcore.xml');

% define exchange reactions IDs
[~, Ex_id] = EditBoundaries(HMRcore, 'Ex_');
IdxExRxns = Ex_id.ID;

% define cooperation reactions IDs
[~, Coop_id] = EditBoundaries(HMRcore, '_COOP');
[~, Biomass_id] = EditBoundaries(HMRcore, 'biomass_synthesis');
IdxCoopRxn = [Coop_id.ID; Biomass_id.ID];


%% Lung Adenocarcinoma 
%% Load single cells and bulk RNAseq data (for different samples
load('data/LUADdataset.mat')

% filter out genes not in HMRcore model
LUAD_filt = extractGeneRow(HMRcore, LUADdataset, 1, DictCORE);

%% Organize in data structures expression profile extracted from LUADdataset data for each cancer samples (H358, LCPT45, LCMBT15). 
% H358 cell line
H358 = makeSCdataset(LUAD_filt.H358_Pooled, LUAD_filt{:,7:56}, LUAD_filt.Properties.VariableNames(7:56), LUAD_filt.HGNC_ID, 10^-4);
H358 = Genes_Sign(H358);
H358 = RepairNegFalse(H358);

% LCPT45 xenograft from primary tumour
LCPT45 = makeSCdataset(LUAD_filt.LCPT45_Pooled, LUAD_filt{:,63:96}, LUAD_filt.Properties.VariableNames(63:96), LUAD_filt.HGNC_ID, 10^-4);
LCPT45 = Genes_Sign(LCPT45);
LCPT45 = RepairNegFalse(LCPT45);

% LCMBT15 zenograft from brain metastasis
LCMBT15 = makeSCdataset(LUAD_filt.LCMBT15_Pooled, LUAD_filt{:,154:202}, LUAD_filt.Properties.VariableNames(154:202), LUAD_filt.HGNC_ID, 10^-4);
LCMBT15 = Genes_Sign(LCMBT15);
LCMBT15 = RepairNegFalse(LCMBT15);

%% Performe Damiani et al integration and analisys
try
    changeCobraSolver('gurobi');
catch
    changeCobraSolver('glpk');
end
%H358
% Compute RAS and integrate them in modelPop
% Modify constraints as in Damiani et al
HMRcore = ScFBAExpSetting(HMRcore, length(H358.CellType));

H358 = single2IntPopModel(H358, HMRcore, IdxExRxns, IdxCoopRxn, 's'); 
fluxIntH358 = optimizeCbModel(H358.modelFVAInt); % optimize model with RAS integrated

% Split single optimize solution in nPop solutions, one for each single
% cells for clustering analisys
[FluxPopH358, RxnPop] = splitScFluxes(H358.modelFVAInt, fluxIntH358, length(H358.CellType));
FluxPopNormH358 = FluxPopH358; % normalize fluxes between 0 and 1.
for i=1:length(RxnPop)
    FluxPopNormH358(i,:) = ((FluxPopH358(i,:) - min(FluxPopH358(i,:)))./(max(FluxPopH358(i,:)) - min(FluxPopH358(i,:))));
end
clst = clustergram(FluxPopNormH358( ~isnan(FluxPopNormH358(:,1)) ,:),'Cluster',2,'Colormap',redgreencmap,'symmetric',false);
addTitle(clst, 'H358 fluxes');

%LCPT45
% Compute RAS and integrate them in modelPop
HMRcore = ScFBAExpSetting(HMRcore, length(LCPT45.CellType));

LCPT45 = single2IntPopModel(LCPT45, HMRcore, IdxExRxns, IdxCoopRxn, 's'); % Compute RAS and integrate them in modelPop
fluxIntLCPT45 = optimizeCbModel(LCPT45.modelFVAInt); % optimize model with RAS integrated

% Split single optimize solution in nPop solutions, one for each single
% cells for clustering analisys
[FluxPopLCPT45, RxnPop] = splitScFluxes(LCPT45.modelFVAInt, fluxIntLCPT45, length(LCPT45.CellType));
FluxPopNormLCPT45 = FluxPopLCPT45; % normalize fluxes between 0 and 1.
for i=1:length(RxnPop)
    FluxPopNormLCPT45(i,:) = ((FluxPopLCPT45(i,:) - min(FluxPopLCPT45(i,:)))./(max(FluxPopLCPT45(i,:)) - min(FluxPopLCPT45(i,:))));
end
clustergram(FluxPopNormLCPT45( ~isnan(FluxPopNormLCPT45(:,1)) ,:),'Cluster',2,'Colormap',redgreencmap,'symmetric',false)
title('LCPT45');

%LCMBT15
% Compute RAS and integrate them in modelPop
HMRcore = ScFBAExpSetting(HMRcore, length(LCMBT15.CellType));

LCMBT15 = single2IntPopModel(LCMBT15, HMRcore, IdxExRxns, IdxCoopRxn, 's'); % Compute RAS and integrate them in modelPop
fluxIntLCMBT15 = optimizeCbModel(LCMBT15.modelFVAInt); % optimize model with RAS integrated

% Split single optimize solution in nPop solutions, one for each single
% cells for clustering analisys
[FluxPopLCMBT15, RxnPop] = splitScFluxes(LCMBT15.modelFVAInt, fluxIntLCMBT15, length(LCMBT15.CellType));
FluxPopNormLCMBT15 = FluxPopLCMBT15; % normalize fluxes between 0 and 1.
for i=1:length(RxnPop)
    FluxPopNormLCMBT15(i,:) = ((FluxPopLCMBT15(i,:) - min(FluxPopLCMBT15(i,:)))./(max(FluxPopLCMBT15(i,:)) - min(FluxPopLCMBT15(i,:))));
end
clustergram(FluxPopNormLCMBT15( ~isnan(FluxPopNormLCMBT15(:,1)) ,:),'Cluster',2,'Colormap',redgreencmap,'symmetric',false)
title('LCMBT15');

%% Breast
%% Load dataset
load('data/BC04dataset.mat')
load('data/BCLN03dataset.mat')

BC04_filt = extractGeneRow(HMRcore, BC04dataset, 1, DictCORE);
BC03LN_filt = extractGeneRow(HMRcore, BC03LNdataset, 1, DictCORE);

% Build scStructures
% BC04 primary breast cancer
BC04 = makeSCdataset(BC04_filt.BC04_Pooled, BC04_filt{:,7:end}, BC04_filt.Properties.VariableNames(7:end), BC04_filt.HGNC_ID, 10^-4);
BC04 = Genes_Sign(BC04);
BC04 = RepairNegFalse(BC04);

% BC03LN lymph node metastasis
BC03LN = makeSCdataset(BC03LN_filt.BC03LN_Pooled, BC03LN_filt{:,7:end}, BC03LN_filt.Properties.VariableNames(7:end), BC03LN_filt.HGNC_ID, 10^-4);
BC03LN = Genes_Sign(BC03LN);
BC03LN = RepairNegFalse(BC03LN);

%% Performe analisys

try
    changeCobraSolver('gurobi');
catch
    changeCobraSolver('glpk');
end

%BC04
% Compute RAS and integrate them in modelPop
HMRcore = ScFBAExpSetting(HMRcore, length(BC04.CellType));

BC04 = single2IntPopModel(BC04, HMRcore, IdxExRxns, IdxCoopRxn, 's'); % Compute RAS and integrate them in modelPop
fluxIntBC04 = optimizeCbModel(BC04.modelFVAInt); % optimize model with RAS integrated

% Split single optimize solution in nPop solutions, one for each single
% cells for clustering analisys
[FluxPopBC04, RxnPop] = splitScFluxes(BC04.modelFVAInt, fluxIntBC04, length(BC04.CellType));
FluxPopNormBC04 = FluxPopBC04; % normalize fluxes between 0 and 1.
for i=1:length(RxnPop)
    FluxPopNormBC04(i,:) = ((FluxPopBC04(i,:) - min(FluxPopBC04(i,:)))./(max(FluxPopBC04(i,:)) - min(FluxPopBC04(i,:))));
end
clustergram(FluxPopNormBC04( ~isnan(FluxPopNormBC04(:,1)) ,:),'Cluster',2,'Colormap',redgreencmap,'symmetric',false)
title('BC04');

%BC03LN
% Compute RAS and integrate them in modelPop
HMRcore = ScFBAExpSetting(HMRcore, length(BC03LN.CellType));

BC03LN = single2IntPopModel(BC03LN, HMRcore, IdxExRxns, IdxCoopRxn, 's'); % Compute RAS and integrate them in modelPop
fluxIntBC03LN = optimizeCbModel(BC03LN.modelFVAInt); % optimize model with RAS integrated

% Split single optimize solution in nPop solutions, one for each single
% cells for clustering analisys
[FluxPopBC03LN, RxnPop] = splitScFluxes(BC03LN.modelFVAInt, fluxIntBC03LN, length(BC03LN.CellType));
FluxPopNormBC03LN = FluxPopBC03LN; % normalize fluxes between 0 and 1.
for i=1:length(RxnPop)
    FluxPopNormBC03LN(i,:) = ((FluxPopBC03LN(i,:) - min(FluxPopBC03LN(i,:)))./(max(FluxPopBC03LN(i,:)) - min(FluxPopBC03LN(i,:))));
end
clustergram(FluxPopNormBC03LN( ~isnan(FluxPopNormBC03LN(:,1)) ,:),'Cluster',2,'Colormap',redgreencmap,'symmetric',false)
title('BC03LN');



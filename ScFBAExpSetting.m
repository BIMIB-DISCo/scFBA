function [model] = ScFBAExpSetting(model,nPops)
% Sets up the boundaries of a given model to respect the choice made in
% scFBA paper. Use it to replicate the resoult.
% 
% USAGE:
%
%   [model] = ScFBAExpSetting(model,nPops)
%
%
% INPUT:
%   model:      metabolic model in COBRA format
%
% OPTIONAL INPUTS:
%   nPops:      if model is a population metabolic model pass number of population.
%               (Default = 1)
%
% OUTPUT:
%   scStruct:   Single Cells dataset with fields added.
%
%
% .. Author:
%       - Chiara Damiani 30/01/2018


if nargin < 2
    nPops = 1;
end

%LOOPS
model = EditBoundaries(model, 'SLC25A5_SLC25A4_SLC25A6', 0,10000, false, false);


%*********************************
%%% vincoli exchange e cooperation reactions
% set boundaries - EX reaction
model = EditBoundaries(model, 'Ex_biomass', 0,10000, false, false);
model = EditBoundaries(model, 'Ex_glucose', -100*nPops,0, false, false);
model = EditBoundaries(model, 'Ex_O2', -100*nPops,0, false, false);
model = EditBoundaries(model, 'Ex_glutamine', -100*nPops,0, false, false);
model = EditBoundaries(model, 'Ex_folate', -100*nPops,0, false, false);
model = EditBoundaries(model, 'Ex_lactateL', 0,0, false, false);
model = EditBoundaries(model, 'Ex_urea', 0,0, false, false);
model = EditBoundaries(model, 'Ex_glutamate', 0,0, false, false);
model = EditBoundaries(model, 'Ex_NH3', 0,0, false, false);
model = EditBoundaries(model, 'Ex_arginines', -100*nPops,0 , false, false);
model = EditBoundaries(model, 'Ex_citratec', 0,0, false, false);
model = EditBoundaries(model, 'Ex_Carnitine', -100*nPops,0, false, false);
model = EditBoundaries(model, 'Ex_Palmitate', 0,0, false, false);
model = EditBoundaries(model, 'Ex_AKG', 0,0, false, false);
model = EditBoundaries(model, 'Ex_proline', 0,0, false, false);
model = EditBoundaries(model, 'Ex_alanine', 0,0, false, false);
model = EditBoundaries(model, 'Ex_aspartate', 0,0, false, false);
model = EditBoundaries(model, 'Ex_asparagine', 0,0, false, false);
model = EditBoundaries(model, 'Ex_serine', 0,0, false, false);
model = EditBoundaries(model, 'Ex_glycine', 0,0, false, false);
model = EditBoundaries(model, 'Ex_cystine', 0,0, false, false);

%set boundaries EX_SC reactions
model = EditBoundaries(model, 'Ex_SC_mercaptopyruvatec', -10000,0, false, false);
model = EditBoundaries(model, 'Ex_SC_H2O', -10000,10000, false, false);
model = EditBoundaries(model, 'Ex_SC_Pi', -10000,10000, false, false);
model = EditBoundaries(model, 'Ex_SC_H', -10000,10000, false, false);

model = EditBoundaries(model, 'DM_glucose', 0,0, false, false);
model = EditBoundaries(model, 'DM_glutamine', 0,0, false, false);
model = EditBoundaries(model, 'DM_folate', 0,0, false, false);
model = EditBoundaries(model, 'DM_lactateL', 0,10000, false, false);
model = EditBoundaries(model, 'DM_urea', 0,10000, false, false);
model = EditBoundaries(model, 'DM_glutamate', 0,10000, false, false);
model = EditBoundaries(model, 'DM_NH3', 0,10000, false, false);
model = EditBoundaries(model, 'DM_arginines', 0,0, false, false);
model = EditBoundaries(model, 'DM_Carnitine', 0,0, false, false);
model = EditBoundaries(model, 'DM_Palmitate', 0,10000, false, false);
model = EditBoundaries(model, 'DM_CO2', 0,10000, false, false);
model = EditBoundaries(model, 'DM_AKG', 0,0, false, false);
model = EditBoundaries(model, 'DM_proline', 0,0, false, false);
model = EditBoundaries(model, 'DM_alanine', 0,0, false, false);
model = EditBoundaries(model, 'DM_aspartate', 0,0, false, false);
model = EditBoundaries(model, 'DM_asparagine', 0,0, false, false);
model = EditBoundaries(model, 'DM_serine', 0,0, false, false);
model = EditBoundaries(model, 'DM_glycine', 0,0, false, false);

% set boundaries - coop reaction (siccome ultimo argomento ? false la
% cambia automatiamente in tutte le cellule)

model = EditBoundaries(model, 'Glucose_DM_COOP', -10000, 0, false, false);
model = EditBoundaries(model, 'Glutamine_DM_COOP', -10000, 0, false, false);
model = EditBoundaries(model, 'Glutamate_DM_COOP', -10000,10000, false, false); %la reazione con consumo di ATP per il demand non viene utilizzata
model = EditBoundaries(model, 'NH3_DM_COOP', -10000,10000, false, false);
model = EditBoundaries(model, 'Palmitate_UP_COOP', 0, 10000, false, false);
model = EditBoundaries(model, 'Palmitate_DM_COOP', 0, 10000, false, false);
model = EditBoundaries(model, 'LactateL_DM_COOP', -10000, 10000, false, false);
model = EditBoundaries(model, 'Oxygen_DM_COOP', -10000,0, false, false);
model = EditBoundaries(model, 'Arginine_DM_COOP', -10000, 0, false, false);
model = EditBoundaries(model, 'Folate_DM_COOP', -10000, 0, false, false);
model = EditBoundaries(model, 'Urea_DM_COOP', 0,0, false, false); %urea
model = EditBoundaries(model, 'Carnitine_DM_COOP', -10000,0, false, false);
model = EditBoundaries(model, 'AKG_DM_COOP', 0,0, false, false);
model = EditBoundaries(model, 'proline_DM_COOP', 0,0, false, false);
model = EditBoundaries(model, 'alanine_DM_COOP', 0,0, false, false);
model = EditBoundaries(model, 'aspartate_DM_COOP', 0,0, false, false);
model = EditBoundaries(model, 'asparagine_DM_COOP', 0,0, false, false);

model = EditBoundaries(model, 'serine_DM_COOP', 0,0, false, false);
model = EditBoundaries(model, 'glycine_DM_COOP', 0,0, false, false);
model = EditBoundaries(model, 'TranspCystineGlu', 0,0, false, false);
model = EditBoundaries(model, 'TranspCystineSer', 0,0, false, false);

%*********************************


%%%reazioni da silenziare
model = EditBoundaries(model, 'Complex1', 0,0, false, false);
model = EditBoundaries(model, 'Complex1ROS', 0,10000, false, false);

%Remove mithocondrial complex I-IV rules
GPR2Remove = {'Complex1','Complex1ROS','Complex2','Complex3','Complex4'};
model.rules(findIdString(model.rxns, GPR2Remove)) = {''};

model.grRules = [];

model = buildRxnGeneMat(model); %generate rxnGeneMat field
model.rev = zeros(length(model.rxns), 1);
model.rev(model.lb < 0 & model.ub > 0) = 1; %generate rev field for reversible reaction

model.metCharges = zeros(length(model.mets), 1);

end

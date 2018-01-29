function [datasetRAS] = getDataScored(model, GeneList, TranscData, ColName, parallel)
% Compute RAS with a dataset as input.
% 
% USAGE:
%
%   [datasetRAS] = getDataScored(model, GeneList, TranscData, ColName, parallel)
%
%
% INPUT:
%   model:          metabolic model in COBRA format
%   GeneList:       vector with genes identifier in the same order of TranscData
%   TranscData:     Dataset genes x sampleas with expression levels. 
%                   Genes must be sorted as GeneList
%
% OPTIONAL INPUTS:
%   ColName:        Identifier for each samples.
%   parallel:       TRUE for use parallel toolbox to speed up the function.
%                   FALSE otherwise. (Default = TRUE)
%
%
% OUTPUTS:
%   datasetRAS:     Dataset reactions x samples with RAS for each reaction.
%                   Nan if there is not a rule associated to that reaction
%                   or all genes have nan as expression level value.

[~, OrdModGen] = ismember(model.genes, GeneList);
GeneList = GeneList(OrdModGen);
TranscData = TranscData(OrdModGen, :);

nRxn = length(model.rxns);
nSam = size(TranscData, 2);

datasetRAS = table();
vettScore = zeros(nRxn,1);
datasetRAS.Reaction = model.rxns;
AllRules = model.rules;
h = waitbar(0,'Data computed: 0%', 'Name','RAS computing');

for k=1:nSam
    TranscData_k = TranscData(:,k);
    if parallel
        parfor i=1:nRxn
            if(strcmp(AllRules(i),''))
                Coeff = NaN;
            else
                rules = AllRules{i};
                Coeff = getScore(rules, TranscData_k);
            end
            vettScore(i) = Coeff;
        end
    else
        for i=1:nRxn
            if(strcmp(AllRules(i),''))
                Coeff = NaN;
            else
                rules = AllRules{i};
                Coeff = getScore(rules, TranscData_k);
            end
            vettScore(i) = Coeff;
        end
    end
    if nargin < 4
        datasetRAS.(strcat('Var',num2str(k))) = vettScore;
    else
        datasetRAS.(ColName{k}) = vettScore;
    end
    waitbar(k/nSam,h,strcat({'Progress: '}, num2str(k/nSam*100, 3), '%'), 'Name','RAS computing');
end
close(h);
end
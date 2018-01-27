function [datasetScored] = getDataScored(model, GeneList, TranscData, ColName, parallel)

% attenzione GeneList deve essere nello stesso ordine di TranscData.
% Suggerisco di passare la colonna del dataset in cui ci sono i geni.

% Sort GeneList and TranscData to match model.genes order
[~, OrdModGen] = ismember(model.genes, GeneList);
GeneList = GeneList(OrdModGen);
TranscData = TranscData(OrdModGen, :);

nRxn = length(model.rxns);
nSam = size(TranscData, 2);

datasetScored = table();
vettScore = zeros(nRxn,1);
datasetScored.Reaction = model.rxns;
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
        datasetScored.(strcat('Var',num2str(k))) = vettScore;
    else
        datasetScored.(ColName{k}) = vettScore;
    end
    waitbar(k/nSam,h,strcat({'Progress: '}, num2str(k/nSam*100, 3), '%'), 'Name','RAS computing');
end
close(h);
%datasetScored = [['Reaction', ItemList]; datasetScored];
end
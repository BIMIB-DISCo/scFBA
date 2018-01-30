function [scStruct] = Genes_Sign(scStruct, CutOff, Hist, ExpXLS, Epsilon, nBins, FitFunction)
% Compute some statistic over the single cell expression profile and
% generate the field geneNameOff. May export the resoult in excel and draw
% an histogram of the distribution for each significative genes.
% 
% USAGE:
%
%   [scStruct] = Genes_Sign(scStruct, CutOff, Hist, ExpXLS, Epsilon, nBins, FitFunction)
%
%
% INPUT:
%   scStruct:       Single Cells dataset in a structure built with makeSCdataset function.
%
% OPTIONAL INPUTS:
%   CutOff:         A gene is defined significative if the mean of its expression 
%                   profile have a ratio less then CutOff between the bulk
%                   expression level. (Default = 0.1)
%   Hist:           If TRUE draw a histogram of each significative genes
%                   expression profile. (Default = FALSE)
%   ExpXLS:         If TRUE export the result in a excel file. (Default = FALSE)
%   Epsilon:        Value added to single cell transcript levels. (Default = minimum of all expression level)
%   nBins:          If Hist = TRUE states the number of bins to draw.
%   FitFunction:    If Hist = TRUE states the function with which fits the distribution of
%                   each significative genes and draw this curve on the histogram. 
%                   (Default = 'kernel')
%
% OUTPUT:
%   scStruct:       Single Cells dataset with fields added.
%
%
% .. Author:
%       - Davide Maspero 30/01/2018

if nargin < 7
    FitFunction = 'kernel';
end

if nargin < 6
    nBins = 30;
end

if nargin < 5
    Epsilon = min(min(scStruct.TPMsc));
end

if nargin < 4
    ExpXLS = false;
end

if nargin < 3
    Hist = false;
end

if nargin < 2
    CutOff = 0.1;
end

FileType = '-dpng';

DataName = inputname(1);

media = mean(scStruct.TPMsc);
scStruct.Stat = abs(media - scStruct.TPMpl)./media;
scStruct.Stat = [media; scStruct.Stat]; %[mean, delta]

offInd = 0;
offMat = zeros(size(scStruct.CellType,1), 1);
signInd = 0;
signMat = zeros(size(scStruct.CellType,1), 1);

for i=1:size(scStruct.Stat, 2)
    if(scStruct.TPMpl(1,i)<=(2*Epsilon) && scStruct.Stat(1, i)<=(2*Epsilon) && min(scStruct.TPMsc(:,i))<=(2*Epsilon))
        offInd = [offInd i];
        offMat = [offMat scStruct.TPMsc(:,i)];
    else
        if(scStruct.Stat(2,i)<=CutOff)
            signInd = [signInd i];
            signMat = [signMat scStruct.TPMsc(:,i)];
        end
    end
end
scStruct.Genes_off = [offInd(:,2:end); offMat(:,2:end)];
scStruct.ZscoreGenSign = [signInd(:,2:end); signMat(:,2:end)];
scStruct.GenNameOff = scStruct.GenName(scStruct.Genes_off(1,:), 1);
Autoscaled = autoscaling(scStruct.ZscoreGenSign(2:end,:));
scStruct.ZscoreGenSign(2:end,:)=Autoscaled;

% Draw an histogram by normalize (Z score) the values.
if Hist
    mkdir(DataName);
    mkdir(strcat(DataName, '\Hist'));
    for i=1:size(Autoscaled, 2)
        histfit2(Autoscaled(:, i), nBins, FitFunction, 'percent');
        axis([-inf,inf,0,0.5]);
        title(strcat(string(scStruct.GenName(scStruct.ZscoreGenSign(1,i))), '--', string(scStruct.Stat(2, scStruct.ZscoreGenSign(1,i)))'));
        FileName = str2mat(strcat(DataName, '\Hist\', string(scStruct.GenName(scStruct.ZscoreGenSign(1,i))')));
        FileName = strrep(FileName, ':', '_');
        print(FileName, FileType);
    end
end
close all
% Esporta in XLS (NB too much columsn may generate an error during the export process.)
if ExpXLS
    mkdir(DataName);
    filename = strcat(DataName, '\', DataName, '_res.xls');
    xlswrite(filename, scStruct.GenName, 'Gene x Cell', 'A2');
    xlswrite(filename, {'Pooled', 'Mean', 'Delta'}, 'Gene x Cell', 'B1');
    xlswrite(filename, [scStruct.TPMpl' scStruct.Stat(1,:)' scStruct.Stat(2,:)'], 'Gene x Cell', 'B2');
    xlswrite(filename, scStruct.CellType', 'Gene x Cell', 'E1');
    xlswrite(filename, scStruct.TPMsc', 'Gene x Cell', 'E2');
    
    if ~isempty(scStruct.Genes_off)
        xlswrite(filename, scStruct.GenName(scStruct.Genes_off(1,:), 1), 'Genes Off', 'A1');
    end
    
    if ~isempty(scStruct.ZscoreGenSign)
        xlswrite(filename, scStruct.GenName(scStruct.ZscoreGenSign(1,:), 1), 'Gene Cutoff Autoscaled', 'A2');
        xlswrite(filename, scStruct.CellType', 'Gene Cutoff Autoscaled', 'B1');
        xlswrite(filename, Autoscaled', 'Gene Cutoff Autoscaled', 'B2');
    end
end

function Matrix = autoscaling(CellXGen)
nEsp=size(CellXGen,1);
nGen=size(CellXGen, 2);

for i=1:nEsp
    for k=1:nGen
        Media=mean(CellXGen(:,k));
        StdDev = std(CellXGen(:,k));
        if StdDev == 0
            Matrix(i,k)=0;
        else
        Matrix(i,k)=(CellXGen(i,k)-Media)/StdDev;
        end
    end
end















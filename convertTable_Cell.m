function [DataOut] = convertTable_Cell(DataIn, header)
% Convert cell to table and vice versa, keeping the header and the first
% column with row labels
%
% USAGE:
%
%   [varOut] = convertTable_Cell(varIn, header)
%
% INPUT:
%   DataIn:          MATLAB table or cell array.
%
% OPTIONAL INPUTS:
%   header:         TRUE if first row contain column labels in cell array 
%                   or if want to keep variable names if a table is passed. 
%                   FALSE if no header needed. (Default = TRUE)
%
% OUTPUTS:
%   varOut:         cell arrays or MATLAB table converted.


if nargin < 2
    header =  true;
end

if istable(DataIn)
    DataOut = NumTable2CellMatrix(DataIn, header);
elseif isnumeric(DataIn)
    DataIn = num2cell(DataIn);
    DataOut = CellMatrix2NumTable(DataIn, false);
elseif iscell(DataIn)
    DataOut = CellMatrix2NumTable(DataIn, header);
else
    warning('Conversion not allowed')
    return
end
end


function [cellMatrixOut] =  NumTable2CellMatrix(tableIn, header)

[r, c] = size(tableIn);
cellMatrixOut = cell(r,c);

for i=1:c
    if isnumeric(tableIn.(i))
        cellMatrixOut(:,i) = num2cell(tableIn.(i));
    else
        cellMatrixOut(:,i) = tableIn.(i);
    end
end
if header
    cellMatrixOut = [tableIn.Properties.VariableNames; cellMatrixOut];
end
end

function [tableOut] =  CellMatrix2NumTable(CellMatrixIn, header)

[~, c] = size(CellMatrixIn);
tableOut = table();

if header
    try
        varName = cellstr(CellMatrixIn(1,:));
    catch
        error('No header found')
        return
    end
    CellMatrixIn(1,:) = []; % prima riga rimossa
end

for i=1:c
    try
        ri = cellstr(CellMatrixIn(:,i));
    catch
        try
            ri = cell2mat(CellMatrixIn(:,i));
        catch
            warning(char(strcat({'conversion error: column '}, num2str(i), {'type not allowed'})));
        end
    end
    if exist('ri', 'var')
        tableOut.(i) = ri;
    end
end
if header
    tableOut.Properties.VariableNames = varName;
end
end







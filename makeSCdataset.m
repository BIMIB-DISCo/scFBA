function scStruct = makeSCdataset(PoolExp, ExpGxC_SC, CellType, Gene_Names, epsilon)

if nargin < 5
    epsilon = 0;
end

if size(CellType, 1)>1
    scStruct.CellType = CellType;
else
    scStruct.CellType = CellType';
end
if size(Gene_Names, 1)>1
    scStruct.GenName = Gene_Names;
else
    scStruct.GenName = Gene_Names';
end
if size(PoolExp, 1)==1
    scStruct.TMPpl = PoolExp;
else
    scStruct.TMPpl = PoolExp';
end
if size(ExpGxC_SC, 1)==size(CellType, 1)
    scStruct.TMPsc = ExpGxC_SC;
else
    scStruct.TMPsc = ExpGxC_SC';
end

scStruct.TMPsc = scStruct.TMPsc + epsilon;

end



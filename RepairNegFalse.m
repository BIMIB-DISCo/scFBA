function [dataset] = RepairNegFalse(dataset, level, epsilon)
epsilon = 10^-4;
count = 0;
%level 1
for i=1:length(dataset.TMPpl)
    if((dataset.TMPpl(i) > epsilon) && (mean(dataset.TMPsc(:,i)) < epsilon*2))
        dataset.TMPsc(:,i) = dataset.TMPpl(i);
        count = count + 1;
    end
end
disp(strcat({'Modificati '}, num2str(count), {' geni'}));

% level 2
if level > 1
    noZero = dataset.TMPsc;
    noZero(noZero <= epsilon) = Inf;
    minNoZero = min(noZero);
    for i=1:size(dataset.TMPsc, 2)
        dataset.TMPsc(noZero(:,i) == Inf, i) =  minNoZero(i);
    end
    dataset.TMPsc(dataset.TMPsc==Inf) = epsilon;
end
end
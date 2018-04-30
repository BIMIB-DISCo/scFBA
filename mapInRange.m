function res = mapInRange(val, minV, maxV, minR, maxR)
%function res = mapInRange(val, minV, maxV, minR, maxR)
%map the value val between [minR maxR].
%minV and maxV the min and max value that valueIn could reach.
%function work also with array or matrix input.

res = minR + (maxR-minR)*(val - minV)./(maxV-minV);

end
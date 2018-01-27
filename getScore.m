function [result] = getScore(rule, transcriptVect)
% geneList and transcriptVect must have the same order
% rule example: (x(index gene:A) & x(index gene:B)) | x(index gene:C)

ruleParsed = parsRule(rule, transcriptVect);

result = resolveRules(ruleParsed);

end

function [rulePars] = parsRule(rule, transcriptVect)
% Replace 'x(posGene)' with corresponding transcript value
i = 0;
rulePars = '';
while i<length(rule)
    i = i + 1;
    if rule(i)=='x'
        posE = strfind(rule(i:end), ')');
        genIndex = str2num(rule(i+2:posE(1)+i-2));
        transcVal = num2str(transcriptVect(genIndex));
        rulePars = [rulePars, transcVal];
        i = posE(1)+i-1;
    else
        rulePars = [rulePars rule(i)];
    end
end

%check rule integrity
if length(strfind(rulePars, ')')) ~= length(strfind(rulePars, '('))
    error(['bracket error on rule: ' rule])
    return
end

end

function [result] = resolveRules(rulePars)

if ~(isempty(strfind(rulePars, '|')) && isempty(strfind(rulePars, '&')) && isempty(strfind(rulePars, '(')) && isempty(strfind(rulePars, ')'))) % controlla che la regola non sia risolta, in tal caso non fa nulla
    
    warning('off', 'all');
    result = NaN;
    pStart = 1;
    i=1;
    while (rulePars(i)~=')' && i < length(rulePars))
        if(rulePars(i)=='(')
            pStart = i;
        end
        i=i+1;
    end
    pEnd = i;
    tmpRule = rulePars(pStart:pEnd);
    tmpRule = strrep(tmpRule, '(', '');
    tmpRule = strrep(tmpRule, ')', '');
    
    subRulesOR = strtrim(split(tmpRule, '|'));
    
    for j=1:size(subRulesOR, 1)
        subRulesAND = strtrim(split(subRulesOR(j), '&'));
        tmpval = nanmin(str2double(subRulesAND));
        if find(~isnan([result tmpval]))
            result = nansum([result tmpval]);
        else
            result = NaN; % se nessun gene esiste restituisce NaN
        end
    end
    %disp(rulePars); % for debug
    rulePars = strrep(rulePars, rulePars(pStart:pEnd), [' ', num2str(result), ' ']);
    result = resolveRules(rulePars); %rilancia funzione ricorsivamente
else
    warning('on', 'all');
    result = str2double(rulePars); %se la regola è risolta, restituisce il valore
end
end
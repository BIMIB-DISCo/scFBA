%USAGE
%findIdString( model.rxns, 'biomass_synthesis')
%findIdString( model.met, 'AKGc')


function [Id_string] = findIdString(where, string)
% funzione che individua l'ID della reazione(singola stringa) o dell reazioni (c)  o del metabolita a partire dal nome

if ~iscell(string)
    string={string};
end

for j=1:length(string)
    Id_string(j)=0;
    
    for i=1:length(where)
        
        
        if strcmp(where(i),string(j))
            
            Id_string(j)=i;
        end
    end
end

end


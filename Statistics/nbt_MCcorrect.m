function Pindex=nbt_MCcorrect(p,Type)
%This function returns the index of MC corrected p-values.

%would be interesting to add more functions..like permutation test

switch Type
    case 'holm'
        %holm correct
        Pc=nbt_holmcorrect(p);
    case 'hochberg'
        %hochberg correct
        Pc = nbt_HochbergCorrect(p);
    case 'bino'
        %bino correct
        if((1-binocdf(length(find(p<0.05)),length(p),0.05))< 0.05)
            Pc = p;
        else
            Pc = [];
        end
    case 'bonfi'
        Pc = p;
        Pc(p > (0.05/length(p))) = nan;
end
if(~isempty(Pc))
    Pindex = nbt_searchvector(p,Pc);
else
    Pindex = [];
end

end
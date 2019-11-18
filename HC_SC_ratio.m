function r = HC_SC_ratio(g)

    hc = 0;
    sc = 0;
    for i=1:length(g.cells)-1
       if g.dead(i)==0
           switch g.populations(i)
               case 2
                   sc = sc + 1;
               case 3
                   hc = hc + 1;
           end
       end
    end

    r = hc/sc;
end
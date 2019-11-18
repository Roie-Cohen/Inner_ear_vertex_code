function g = initiateDeltaSignal(g)

    for i=1:length(g.cells)-1
        switch g.populations(i)
            case 2
                g.signal(i) = 0;
            case 3
                g.delta(i) = 18 + 4*rand;
        end
    end
end
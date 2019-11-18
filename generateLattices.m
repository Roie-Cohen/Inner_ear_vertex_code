%% generating lattice
dir = 'lattices';
mkdir(dir)
nlat = 1; %% number of lattices
regl = 'ri';
lt = 0; 
for reg = lt:lt %% if reg==1, regular lattice (6 neighbor/cell) or not. 
    allper = [];
    bfname = [dir '/' regl(reg+1) 'lat12x12_'];
    for nl = 1:nlat,
        close all;
        g = hexIrregLattice(12,12,reg);
%         g = ImportedLattice(reg);
        fname = [bfname num2str(nl)]
        print(gcf,'-depsc2',[fname '.eps']);    % Encapsulated Level 2 Color PostScript
        if lt == 3,
            [weights perim area] = getConnectivity(g);
            save(fname,'g','weights','perim','area');
            allper = [allper perim];
        end
    end
    figure, hist(allper);
    xlabel('perimeter');
    print(gcf,'-depsc2',[bfname 'hist.eps']);    % Encapsulated Level 2 Color PostScript
end


%% reloading lattices
if(0),
   bfname = [dir 'lattices/lat12x12_'];
   for nl = [2:20],
        close all;
         fname = [bfname num2str(nl)];
         load (fname);
         LatticePresentation(g,0);
         print(gcf,'-depsc2',[fname '.eps']);
   end
end
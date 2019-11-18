% load some disordered lattice
load('lattices/diff_lat2');

cc = 1:length(g.cells)-1;
c = length(cc);
            
% removing double bonds apearing in the bonds array
bonds = g.bonds(g.bonds(:,1) ~= 0, :);
for i=1:length(bonds)
    mult = find( sum(bonds(:,1:2) == flip(bonds(i,1:2),2), 2) == 2 ); % finding the multiplicity index
    if ~isempty(mult)
        bonds(i,:) = 0;
    end
end
bonds = bonds(bonds(:,1) ~= 0, :); %at this point there are no multiple bonds in the array
nb = length(bonds(:,1));

% we assume that the area in linearly proportional to the pressure
apf = 1;    % the area-pressure factor
pres = apf*g.areas';

% finding the relevent vertices
rverts = unique([bonds(:,1); bonds(:,2)]);
rv_pos = getRelativePosition(g, rverts);

np = nb + length(pres);     % number of parameters
nc = 2*length(rverts);      % number of conditions

% building the M matrix
M = zeros(nc+1,np);
for i=1:length(rverts)
   vi = rverts(i);   % vertex index
   [Trow, Tcol] = find( bonds(:,1:2) == vi );
   Tcol = 3 - Tcol; % changing to the neighbor index (1->2, 2->1)
   nbrs = zeros(length(Trow),1); % the vertex neighbors
   for n=1:length(nbrs)
       nbrs(n) = bonds(Trow(n), Tcol(n));
   end
   
   for n=1:length(nbrs)
       vj = nbrs(n);    % neighbor index
       vij = getRelativePosition(g, vj) - getRelativePosition(g, vi);
       rij = norm(vij);
       
       bcells = [bonds(Trow,3), bonds(Trow,4)]; % bordering cells by the i-j bond
       % if Tcol(n) == 2 then the bordering cells are given in a clockwise
       % order relative to the spatial vector going from i to j (vij).
       % if not we'll flip it to get clockwise order.
       if Tcol(n) == 1
           bcells = flip(bcells);
       end
       
       % x component of the force exerted by j
       M(2*i-1, Trow(n)) = vij(1)/rij;          % the factor before the tention coefficient
       M(2*i-1, nb + bcells(1)) = -0.5*vij(2);  % the factor before the first cell pressure coefficient
       M(2*i-1, nb + bcells(2)) = 0.5*vij(2);   % the factor before the second cell pressure coefficient
       % y component of the force exerted by j
       M(2*i, Trow(n)) = vij(2)/rij;            % the factor before the tention coefficient
       M(2*i, nb + bcells(1)) = 0.5*vij(1);     % the factor before the first cell pressure coefficient
       M(2*i, nb + bcells(2)) = -0.5*vij(1);    % the factor before the second cell pressure coefficient
   end
end

M(nc+1,:) = 1;
% additional constraints
add_con = 1;
if add_con
    % constant pressure
    Add = zeros(c,np);
    for i=1:c
        Add(i,nb+i) = 1;
    end
    M = [M; Add];
end

Minv = pinv(M);
C = zeros(nc+1,1);  % condition vector
C(nc+1) = nb;       % scale factor
if add_con
   C = [C; zeros(c,1)]; 
end
psi = Minv*C;


% creating the real tension vector for comparison
if(isfield(g,'populations'))
    rt = zeros(nb,1);
    for i=1:nb
        pp = [g.populations(bonds(i,3)) g.populations(bonds(i,4))];
        if sum(pp == [1 2])==2 || sum(pp == [1 3])==2 || sum(pp == [2 1])==2 || sum(pp == [3 1])==2
            rt(i) = 3;
        else
            if sum(pp == [3 2])==2 || sum(pp == [2 3])==2
                rt(i) = 1;
            else
                if (sum(pp == [3 3])==2) || sum(pp == [2 2])==2
                    rt(i) = 1.5;
                else
                    rt(i) = 1;
                end
            end
        end
    end
end
rt = nb/sum(rt).*rt;    % renormalizing so the sum of rt is nb


% ploting real tension vs infered tension (in absolute value)
it = psi(1:length(rt));
scatter(rt,abs(it))
xlabel('Real tension');
ylabel('Infered tension');

% comparing the average of the indered values to the real tension
ut = unique(rt);
for i=1:length(ut)
    avg = mean(it(rt==ut(i)));
    s = ['Real tension: ', num2str(ut(i)), '    Infered average:', num2str(avg)];
    disp(s);
end

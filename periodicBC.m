function ng= periodicBC(g,nrow,ncol)
ng = g;
bf = ncol*(nrow-1)
%% four corners
%%bottom left
idx = find(ng.bonds(ng.cells{2},4)==0);
 dd =diff(idx);
 st = find(dd~=1);
 if(length(st)>0),
     idx = [idx(st+1:length(idx)); idx(1:st)];
 end
 ng.bonds(ng.cells{2}(idx),4) = [bf+1;ncol*nrow;ncol;2*ncol];
 %% bottom right
 idx = find(ng.bonds(ng.cells{ncol+1},4)==0);
 dd =diff(idx);
 st = find(dd~=1);
 if(length(st)>0),
     idx = [idx(st+1:length(idx)); idx(1:st)]
 end
 ng.bonds(ng.cells{ncol+1}(idx),4) = [1;ncol*nrow;ncol*nrow-1];
 
 %%top left
idx = find(ng.bonds(ng.cells{bf+2},4)==0);
 dd =diff(idx);
 st = find(dd~=1);
 if(length(st)>0),
     idx = [idx(st+1:length(idx)); idx(1:st)];
 end
 ng.bonds(ng.cells{bf+2}(idx),4) = [ncol*nrow;1;2];
 
 %%top right
idx = find(ng.bonds(ng.cells{ncol*nrow+1},4)==0);
 dd =diff(idx);
 st = find(dd~=1);
 if(length(st)>0),
     idx = [idx(st+1:length(idx)); idx(1:st)];
 end
 ng.bonds(ng.cells{ncol*nrow+1}(idx),4) = [ncol;1;bf+1;bf+1-ncol];
 
 
%% first row
for i=2:ncol-1,
   % i
idx = find(ng.bonds(ng.cells{i+1},4)==0);
if(idx(2)-idx(1)>1)
    idx = [idx(2),idx(1)];
end
ng.bonds(ng.cells{i+1}(idx),4) = [bf+i;bf+i-1]; 
end

%% last row
for i=2:ncol-1,
idx = find(ng.bonds(ng.cells{bf+i+1},4)==0);
if(idx(2)-idx(1)>1)
    idx = [idx(2),idx(1)];
end
ng.bonds(ng.cells{bf+i+1}(idx),4) = [i;i+1]; 
%% merging vertices
% for j=1:2,
%     bb = find(ng.bonds(ng.cells{i+j},4)== bf+i); %% i+j is a shortcut
%     ng.bonds(ng.cells{bf+i+1}(idx(j)),1:2) = ng.bonds(ng.cells{i+j}(bb),[2 1]);
% end
end

%%first column
for i=ncol+1:ncol:bf,
   idx = find(ng.bonds(ng.cells{i+1},4)==0);
   if(length(idx)==1),
       ng.bonds(ng.cells{i+1}(idx),4)=i+ncol-1;
   else
       dd =diff(idx);
%        if(dd(1)>1),
%            idx = [idx(2:3); idx(1)];
%        end
%        if (dd(2)>1),
%             idx = [idx(3); idx(1:2)];
%        end
    st = find(dd~=1);
 if(length(st)>0),
     idx = [idx(st+1:length(idx)); idx(1:st)];
 end
      ng.bonds(ng.cells{i+1}(idx),4)=[i-1;i+ncol-1;i+2*ncol-1]; 
   end
end
%% last column
for i=2*ncol:ncol:bf,
   idx = find(ng.bonds(ng.cells{i+1},4)==0);
   if(length(idx)==1),
       ng.bonds(ng.cells{i+1}(idx),4)=i-ncol+1;
      % bb = find(ng.bonds(ng.cells{i+1},4)==i+ncol-1;);
       %% merging vertex
      %% ng.bonds(ng.cells{i+ncol}(idx),[1:2])=  ng.bonds(ng.cells{i+1}(bb),[2 1]);
   else
       dd =diff(idx);
       if(dd(1)>1),
           idx = [idx(2:3); idx(1)];
       end
       if (dd(2)>1),
            idx = [idx(3); idx(1:2)];
       end
      ng.bonds(ng.cells{i+1}(idx),4)=[i+1;i-ncol+1;i-2*ncol+1]; 
   end
end
% ng.xboundary([1:ncol:bf+1],1)=1;
% ng.yboundary([1:ncol],1)=1;
% ng.xboundary([ncol:ncol:ncol*nrow],2)=1;
% ng.yboundary([bf+1:ncol*nrow],2)=1;
% 
% 
% ng.xboundary([2:ncol:bf+2],1)=1;
% ng.yboundary([ncol+1:2*ncol],1)=1;
% ng.xboundary([ncol-1:ncol:ncol*nrow-1],2)=1;
% ng.yboundary([bf+1-ncol:bf],2)=1;
ng.bc=1;
%%% merge vertices
if(1)

%% g.bonds(find(g.bonds(:,4)==0),:)
for i=ncol:-1:1,
%for i=[1],
    for j = length(ng.cells{bf+i+1}):-1:1,
        neb = ng.bonds(ng.cells{bf+i+1}(j),4);
        if(neb>0),
        revb = find(ng.bonds(ng.cells{neb+1},4)==bf+i);
      %  disp(['neb: ' num2str(bf+i) ' ' num2str(neb)])
      %  disp([num2str(ng.bonds(ng.cells{bf+i+1}(j),1:2)) ' or ' num2str(ng.bonds(ng.cells{neb+1}(revb),[2 1]))]);
        aa = ng.bonds(:,1:2);
        aa(aa== ng.bonds(ng.cells{bf+i+1}(j),1))= ng.bonds(ng.cells{neb+1}(revb),2);
        aa(aa== ng.bonds(ng.cells{bf+i+1}(j),2))= ng.bonds(ng.cells{neb+1}(revb),1);
      aa;
        ng.bonds(:,1:2)=aa;
         % ng.bonds(ng.cells{bf+i+1}(j),1:2) = ng.bonds(ng.cells{neb+1}(revb),[2 1]);
        end
        end
end

% last column
for i=ncol:ncol:nrow*ncol,
     for j = 1:length(ng.cells{i+1}),
     neb = ng.bonds(ng.cells{i+1}(j),4);
      revb = find(ng.bonds(ng.cells{neb+1},4)==i);
      if(ng.bonds(ng.cells{i+1}(j),1:2) ~= ng.bonds(ng.cells{neb+1}(revb),[2 1])),
          aa = ng.bonds(:,1:2);
        aa(aa== ng.bonds(ng.cells{i+1}(j),1)) = ng.bonds(ng.cells{neb+1}(revb),2);
         aa(aa== ng.bonds(ng.cells{i+1}(j),2)) = ng.bonds(ng.cells{neb+1}(revb),1);
         ng.bonds(:,1:2)=aa;
      end
     end
end

%% removing unused boundaries
torem = find(ng.bonds(:,3)==0);
ng = remove_bond(ng,torem);

%% rescaling
scaley = 2*pi/(nrow*1.5);
scalex = 2*pi/(ncol*sqrt(3));
ng.verts(:,1) = (g.verts(:,1)-sqrt(3)*1.5)*scalex - pi;
ng.verts(:,2) = (g.verts(:,2)-2)*scaley - pi;

end
function cv = findConnectingVert(skel, vidx, vert, branch)

    brow = vidx(vert,1);
    bcol = vidx(vert,2);
    probe = zeros(3,3);
    probe(branch(1), branch(2)) = 1;
    
    while (1==1)
        [r, c] = find(probe);
        switch numel(r)
            case 0
                cv = 0;
                break;
            case 1
                brow = brow + r - 2;
                bcol = bcol + c - 2;
                probe = [ [skel(brow-1,bcol-1), skel(brow-1,bcol), skel(brow-1,bcol+1)];
                        [  skel(brow  ,bcol-1), 0                , skel(brow  ,bcol+1)];
                        [  skel(brow+1,bcol-1), skel(brow+1,bcol), skel(brow+1,bcol+1)] ];
                probe(4-r, 4-c) = 0;
            otherwise
                cv = find(sum(ismember(vidx,[brow bcol]),2)==2);
                break;
        end
        
    end
end
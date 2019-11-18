function LatticePresentation(g,num, fig_num)
%% Draws the lattice g. If num is true, display the cell numbers
if nargin==1,
    num=0;
end

if nargin==3
    clf(fig_num);
end
if(isfield(g,'delta'))
    isdelta = 0;
    maxd = max(g.delta);
    maxs = max(g.signal);
else
    isdelta = 0;
end
cells = [1:length(g.cells)-1];
edge_color = 'k';
for i = cells
    if(length(g.cells{i+1})>2),
        if(g.dead(i)==0),
            verts = g.bonds(g.cells{i+1}(:),1);
            if(length(verts)>2),
                v = getRelativePosition(g,verts,i);
                if(isfield(g,'scale')),
                    v = v*g.scale;
                end
                if(isfield(g,'populations')),
                    if (isdelta)
                        if (g.delta(i) == 0) || (g.signal(i) == 0)
                            col = [0 0 0];
                        else
                            col = [g.delta(i)/maxd, g.signal(i)/maxs,0];
                        end
                        patch(v(:,1),v(:,2),col,'EdgeColor','k');
                    else
                        switch g.populations(i)
                            case 0
                                patch(v(:,1),v(:,2),'w','EdgeColor','r','LineWidth',1.5);
                            case 1
                                patch('XData',v(:,1), 'YData' ,v(:,2),'FaceColor',[1 0.8 0.8],'EdgeColor',edge_color,'LineWidth',1.5);
                            case 2
                                patch('XData',v(:,1), 'YData' ,v(:,2),'FaceColor',[0.9 0.9 0.9],'EdgeColor',edge_color,'LineWidth',1.5);
                            case 3
                                patch('XData',v(:,1), 'YData', v(:,2),'FaceColor',[0 0.75 0],'EdgeColor',edge_color,'LineWidth',1.5);
                            case 4
                               patch('XData',v(:,1), 'YData', v(:,2),'FaceColor',[0.5 0.5 0.5],'EdgeColor',edge_color,'LineWidth',1.5);
                            case 5
                               patch('XData',v(:,1), 'YData', v(:,2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',edge_color,'LineWidth',1.5);
                        end
                    end
                else
                    patch(v(:,1),v(:,2),'w','EdgeColor','k');
                end
                hold on;
                switch num
                    case 1
                        text(mean(v(:,1)),mean(v(:,2)),num2str(i),'HorizontalAlignment','center','color','k');
                    case 2
                        as = 0.5*sqrt(cellarea(g,i)); %arrow size
                        ax = as*cos(g.polarity(i));
                        ay = as*sin(g.polarity(i));
                        p1 = [mean(v(:,1))-ax/2, mean(v(:,2))-ay/2];
                        p2 = [mean(v(:,1))+ax/2, mean(v(:,2))+ay/2];
                        dp = p2-p1;
                        quiver(p1(1),p1(2),dp(1),dp(2),'b','MaxHeadSize',0.5,'LineWidth',1);
                end
            end
        end
    end
end

axis equal
set(gca,'visible','off');
hold off;
end

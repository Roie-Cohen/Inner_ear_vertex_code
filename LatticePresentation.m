function LatticePresentation(g,num, fig_num)
%% Draws the lattice g. If num is true, display the cell numbers
if nargin==1
    num=0;
end

if nargin==3
    clf(fig_num);
end
cells = [1:length(g.cells)-1];
edge_color = 'k';
for i = cells
    if(length(g.cells{i+1})>2)
        if(g.dead(i)==0)
            verts = g.bonds(g.cells{i+1}(:),1);
            if(length(verts)>2)
                v = getRelativePosition(g,verts,i);
                if(isfield(g,'scale'))
                    v = v*g.scale;
                end
                if(isfield(g,'populations'))
                    switch g.populations(i)
                        case 0
                            patch(v(:,1),v(:,2),'w','EdgeColor',edge_color,'LineWidth',1.5);
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
                else
                    patch(v(:,1),v(:,2),'w','EdgeColor',edge_color);
                end
                hold on;
                if num == 1
                    text(mean(v(:,1)),mean(v(:,2)),num2str(i),'HorizontalAlignment','center','color','k');
                end
            end
        end
    end
end

axis equal
set(gca,'visible','off');
hold off;
end

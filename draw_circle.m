function h = draw_circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
% h = plot(xunit, yunit, 'k', 'linewidth', 2);
h = patch('XData',xunit, 'YData', yunit,'FaceColor',[0 0.75 0],'EdgeColor','r','LineWidth',1.5);
%% figure settings
function Plotter(titlename,xname,yname,xlim,ylim,legendArrary,flag,FontSize)
if nargin<8
    FontSize=12;
end
xlabel(xname);
ylabel(yname);
title(titlename);
grid on;
 set(gca, 'GridLineStyle', '--'); % 设置网格线为虚线
box on;


width = 580;
height = width * 0.7; % 假设高度与宽度比例为 0.7
% axis square;
set(gca, 'FontName', 'Arial', 'FontSize', FontSize);
set(gcf,'Position', [200, 200, width, height]);
set(gca, 'LineWidth', 1.25);
set(gca, 'XLim',xlim,'YLim',ylim);
% ylim tight;
if flag.LegendON_OFF
    lgd = legend(legendArrary{:},'Location','best');
    set(lgd,'FontName','Arial','FontSize',FontSize);
    if flag.Legendflage
        set(lgd, 'Color', 'none'); % 设置图例框的颜色为'none'
        set(lgd, 'Box', 'off');
    end
end
end
% lgd = legend('Experimental','Simulated','Experimental R_p','Location','best');
% set(lgd, 'Color', 'none'); % 设置图例框的颜色为'none'
% set(lgd, 'Box', 'off');
% yticks([-5, -4, -3, -2,-1, 0, 1,2,3,4,5]);

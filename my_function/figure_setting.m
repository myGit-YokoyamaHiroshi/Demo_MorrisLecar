function figure_setting(X, Y, fig)
    xSize = X - 2;
    ySize = Y - 2;

    set(fig, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)
    movegui(fig, 'center');

    % figure size printed on paper
%     set(fig, 'PaperPositionMode','auto')
    
    set(fig, 'PaperUnits','centimeters')
    set(fig, 'PaperSize',[X Y]/2)
    set(fig, 'PaperPosition',[1 1 xSize ySize]/2)
    set(fig, 'PaperOrientation','portrait')
end
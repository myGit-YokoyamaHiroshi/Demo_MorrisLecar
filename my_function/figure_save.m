function figure_save(fig, fname)
    set(gca, 'FontName', 'Arial')
    print(fig, '-depsc2',[fname,'.eps'], '-r300');
    print(fig, '-painters', '-dpdf', [fname,'.pdf'], '-r300');
    print(fig, '-dpng',[fname,'.png'], '-r300');
end
% plot export function
function exportPlot(fileName, fig)
    print(fig, sprintf('%s.pdf', fileName), '-dpdfcairo');
    system(sprintf('mv %s.pdf plots/%s.pdf', fileName, fileName));
end
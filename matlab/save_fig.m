function save_fig(f, figdir, name)
% SAVE_FIG  Save figure f as figdir/name.png.
if ~exist(figdir, 'dir')
    mkdir(figdir);
end
print(f, fullfile(figdir, [name '.png']), '-dpng', '-r150');
end

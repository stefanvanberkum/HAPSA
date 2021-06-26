%% Clear all data.
clear
close all

%% Load the results.
path = "C:/Users/Stefan van Berkum/git/HAPSA/HAPSA/Results/";
filenames = ["APSA/CSV/30_240", "VIS/CSV/30_240_50.0", "VIS/CSV/30_240_100.0", "VIS/CSV/30_240_150.0", "VIS/CSV/30_240_200.0", "HAPSA/CSV/30_240_100.0_5.0E-4"];
outpath = "C:/Users/Stefan van Berkum/Google Drive/Studie/Thesis/Figures/Colormaps/";
outfiles = ["APSA_30_240", "VIS_30_240_50" , "VIS_30_240_100", "VIS_30_240_150", "VIS_30_240_200", "HAPSA_30_240_100_00005"];

%% Generate colormaps.
x = [1, 3];
y = [1, 6];

for i=1:length(filenames)
    matrix = readmatrix(path + filenames(i) + "_shelf.csv");
    plot = imagesc(x, y, flipud(matrix));
    colormap(flipud(summer));
    colorbar;
    caxis([40, 65]);
    set(gca, 'YDir', 'normal');
    set(gca, 'xtick', [1, 2, 3]);
    filename = outpath + outfiles(i);
    saveas(plot, filename, 'epsc')
end

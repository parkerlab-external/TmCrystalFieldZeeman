pumpprobe_files = [
    ".\vspectrum2024-09-30-1735_processed.mat";
    ".\vspectrum2024-10-02-2054_processed.mat";
    ];

PPfig = figure("Name","PumpProbe", "Units","inches", "Position", [0 0 5.25 2.5]);

axis1 = axes(PPfig);
axes(axis1);
hold on;

%First pumpprobe site - second index for better resolution on the probe
load(pumpprobe_files(1));
plot(axis1, scan_wvs, probe+380, "black");

site1col = "#db4939";

nonlin = pump + probe - pumpprobe;
plot(axis1, scan_wvs, 1.5*(nonlin - mean(nonlin)) + 220, "Color", site1col);
text(1140.9, 220, "$\times 1.5$", 'Interpreter','latex');
xline(pump_wv, "Alpha",0.3, "LineWidth",2, "Color", site1col,...
    "Label",sprintf("%.3f", pump_wv),"FontName","Times New Roman", "LabelHorizontalAlignment","left");

%Second pumpprobe site
load(pumpprobe_files(2));
site2col = "#3c70c9";

nonlin = pump + probe - pumpprobe;
plot(axis1, scan_wvs, 3*(nonlin - mean(nonlin)) + 12, "Color", site2col);
text(1140.9, 30,"$\times 3$","interpreter","latex");
xline(axis1, pump_wv, "Alpha",0.3, "LineWidth",2, "Color",site2col,...
    "Label",sprintf("%.3f", pump_wv),"FontName","Times New Roman", "LabelHorizontalAlignment","left");

axis1.XLim = [1139 1141];
% axis1.Color = "none";

xlabel("Probe Wavelength (nm)")
ylabel(sprintf("Fluorescence Amplitude\n(counts/ms)"));
set(axis1, 'FontName','Times New Roman');

axis2 = axes;
plot(axis2, scan_wvs, probe+380);
axis2.XAxisLocation = "top";
axis2.YAxisLocation = "right";
axis2.Color = "none";
axis2.XLim = 1e7./[axis1.XLim(2) axis1.XLim(1)];


axis2.XLabel.String = "Energy (cm^{-1})";
axis2.Position = axis1.Position;
axis2.Box = "off";

axis1.YLim = [-100 800];
axis2.YTick = axis1.YTick;
axis2.YTickLabel = [];
axis2.YLim = axis1.YLim;

%Janky part: steal the axis from the auto axis creation
[autopositions, sortind] = sort(1e7./axis2.XTick);
autolabels = axis2.XTickLabel(sortind);
axis2.XLim = axis1.XLim;
axis2.XTick = autopositions;
axis2.XTickLabel = autolabels;

set(axis2, 'FontName', 'Times New Roman');

set([axis1, axis2], "Position", [0.12, 0.14, 0.85, 0.71])
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 10);

axis1.XRuler.TickLabelGapOffset = 0;
axis2.XRuler.TickLabelGapOffset = 0;

load(fullfile("SiteILevelFitted.mat"))
line_positions = unique(eE) - unique(gE)';
line_positions = 1e7./line_positions(:);
num_lines = numel(line_positions);
for i = 1:num_lines
    xpos = axis1.Position(1) + axis1.Position(3)*(line_positions(i)-axis1.XLim(1))/(axis1.XLim(2) - axis1.XLim(1));
    if xpos < 1
        annotation(PPfig, 'arrow', [xpos xpos], [0.34 0.37], "Color", site1col, "HeadLength",3, "HeadWidth",3);
    end
end

load(fullfile("SiteIILevelFitted.mat"));
line_positions = unique(eE) - unique(gE)';
line_positions = 1e7./line_positions(:);
num_lines = numel(line_positions);
for i = 1:num_lines
    xpos = axis1.Position(1) + axis1.Position(3)*(line_positions(i)-axis1.XLim(1))/(axis1.XLim(2) - axis1.XLim(1));
    if xpos < 1
        annotation(PPfig, 'arrow', [xpos xpos], [0.17 0.2], "Color", site2col, "HeadLength",3, "HeadWidth",3);
    end
end


figure("Name", "Levels","Units","inches", "Position",[5.25,0,1.8,2.5])
myax = axes();
load(fullfile("SiteILevelFitted.mat"));
gE = unique(gE);
eE = unique(eE) - 8771.23;
for i= 1:numel(gE)
    line([0 1], [gE(i) gE(i)], "Color", site1col);
end
for i= 1:numel(eE)
    line([3 4], [eE(i) eE(i)], "Color", site1col);
end
load(fullfile("SiteIILevelFitted.mat"));
gE = unique(gE);
eE = unique(eE) - 8771.23;
for i= 1:numel(gE)
    line([1 2], [gE(i) gE(i)], "Color", site2col);
end
for i= 1:numel(eE)
    line([4 5], [eE(i) eE(i)], "Color", site2col);
end
% myax.Color = "none";
% myax.XTick = [];
% myax.YTick = [];
axis off





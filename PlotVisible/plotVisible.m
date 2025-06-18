colors = ["cyan"
    "green"
    "blue"
    "red"];

figure("Name","PlotVisible");
clf;
vis_ax = axes();
hold on;

load("vspectrum2024-10-09-1341.mat");
pump_wvs = [];
for i = [1 4]
    mean_data = squeeze(mean(data_out.raw_data(i,:,:,:),3));
    mean_data = mean(mean_data(:,1:75),2);
    plot(vis_ax, data_out.emission_range, mean_data, "Color", colors(i),"LineWidth",1);
    pump_wvs = [pump_wvs shift_wavelengths(data_out.wv_samples(i))];
end
hold off
legend(string(pump_wvs))

xlim(vis_ax, [data_out.emission_range(1) 800]);
set(vis_ax, "FontName", 'Times New Roman');
vis_ax.XLabel.String = "Wavelength (nm)";
vis_ax.YLabel.String = "Luminescence (counts/ms)";

axis2 = axes();
axis2.XAxisLocation = "top";
axis2.YAxisLocation = "right";
axis2.YTick = vis_ax.YTick;
axis2.YTickLabel = [];
axis2.YLim = vis_ax.YLim;
axis2.XLim = sort(1e7./vis_ax.XLim);
axis2.XDir = "reverse";
axis2.XLabel.String = "Energy (cm^{-1})";
axis2.Position = vis_ax.Position;
axis2.Box = "off";
axis2.Color = "none";
set(axis2, "FontName", 'Times New Roman');

title(axis2, "Emission");
set([vis_ax, axis2],"Position", [0.13 0.11 0.775 0.75]);



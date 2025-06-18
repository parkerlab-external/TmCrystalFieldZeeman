load("PolarizedSpectra.mat");

PSfig = figure("Name","Polarized Spectra","Units", "inches","Position", [0.5 0.5 3.4 3.3]);
main_ax = axes();
main_ax.Box = "on";

%Sort
[scan_wvs1, sortind] = sort(scan_wvs1);
scan1 = scan1(sortind);
[scan_wvs2, sortind] = sort(scan_wvs2);
scan2 = scan2(sortind);
[scan_wvs3, sortind] = sort(scan_wvs3);
scan3 = scan3(sortind);

plot(main_ax, scan_wvs1, scan1, "Color", "#a836c7");
hold on;
plot(main_ax, scan_wvs2, scan2, "black");
plot(main_ax, scan_wvs3, scan3, "black",LineStyle='-.', LineWidth=1);

xlimits = [1140.01 1140.0165];

%integrated strengths
m = scan_wvs1 <= xlimits(2) & scan_wvs1 >= xlimits(1);
IS1 = trapz(scan_wvs1(m), scan1(m) - min(scan1(m)));
m = scan_wvs2 <= xlimits(2) & scan_wvs2 >= xlimits(1);
IS2 = trapz(scan_wvs2(m), scan2(m) - min(scan2(m)));
m = scan_wvs3 <= xlimits(2) & scan_wvs3 >= xlimits(1);
IS3 = trapz(scan_wvs3(m), scan3(m) - min(scan3(m)));
IS = mean([IS1 IS2 IS3]);

linewidth = 0.001;
Bfield = 0.01;
scale = 1/2;

xaxisinput = 1e7./(linspace(xlimits(1), xlimits(2), 300));
Ninput = 100;

%M1 only
M1color = "#db4939";
[xaxis, predicted] = FullInteraction(Bfield, linewidth, "parallel1", "fittedPQ.mat", 1, 0, "N",Ninput, "xaxis", xaxisinput);
[xaxis, sortind] = sort(1e7./xaxis);
predicted = predicted(sortind);
factor = IS/trapz(xaxis, predicted)*scale;
plot(main_ax, xaxis, factor*predicted - 15, Color=M1color);

[xaxis, predicted] = FullInteraction(Bfield, linewidth, "transverse1", "fittedPQ.mat", 1, 0, "N",Ninput,"xaxis", xaxisinput);
[xaxis, sortind] = sort(1e7./xaxis);
predicted = predicted(sortind);
plot(main_ax, xaxis, factor*predicted - 15, Color=M1color, LineStyle='-.', LineWidth=1);


%PQOnly
PQcolor = "#3c70c9";
[xaxis, predicted] = FullInteraction(Bfield, linewidth, "parallel1", "fittedPQ.mat", 0, 1, "N",Ninput,"xaxis", xaxisinput);
[xaxis, sortind] = sort(1e7./xaxis);
predicted = predicted(sortind);
factor = IS/trapz(xaxis, predicted)*scale;
plot(main_ax, xaxis, factor*predicted - 30, Color=PQcolor);

[xaxis, predicted] = FullInteraction(Bfield, linewidth, "transverse1", "fittedPQ.mat", 0, 1,"N",Ninput,"xaxis", xaxisinput);
[xaxis, sortind] = sort(1e7./xaxis);
predicted = predicted(sortind);
plot(main_ax, xaxis, factor*predicted - 30, Color=PQcolor, LineStyle='-.', LineWidth=1);

%FDOnly
FDcolor = "#36c781";
[xaxis, predicted] = FullInteraction(Bfield, linewidth, "parallel1", "fittedFD.mat", 0, 1,"N",Ninput,"xaxis", xaxisinput);
[xaxis, sortind] = sort(1e7./xaxis);
predicted = predicted(sortind);
factor = IS/trapz(xaxis, predicted)*scale;
plot(main_ax, xaxis, factor*predicted - 45, Color=FDcolor);

[xaxis, predicted] = FullInteraction(Bfield, linewidth, "transverse1", "fittedFD.mat", 0, 1,"N",Ninput,"xaxis", xaxisinput);
[xaxis, sortind] = sort(1e7./xaxis);
predicted = predicted(sortind);
plot(main_ax, xaxis, factor*predicted - 45, Color=FDcolor, LineStyle='-.', LineWidth=1);

annotation(PSfig, "textbox", [0.5 0.5 0.1 0.1],"String","Data", 'FitBoxToText', "on", "FontSize",10, "FontName","Times New Roman");
annotation(PSfig, "textbox", [0.5 0.4 0.1 0.1],"String","M1", 'FitBoxToText', "on", "FontSize",10, "FontName","Times New Roman");
annotation(PSfig, "textbox", [0.5 0.3 0.1 0.1],"String","PQ", 'FitBoxToText', "on", "FontSize",10, "FontName","Times New Roman");
annotation(PSfig, "textbox", [0.5 0.2 0.1 0.1],"String","FD", 'FitBoxToText', "on", "FontSize",10, "FontName","Times New Roman");

legend(["" "E || B" "E \perp B"])

set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','FontName'),'FontName',"Times New Roman");

main_ax.XLim = [1140.01 1140.016];
main_ax.XLabel.String = "Wavelength (nm)";
main_ax.YLim = [-50 40];
main_ax.YLabel.String = "Fluorescence Amplitude (counts/ms)";





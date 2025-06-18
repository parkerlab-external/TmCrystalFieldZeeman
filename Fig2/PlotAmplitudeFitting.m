%Load measured data
load(fullfile("..","Saves","rough_position_scan"));

scan_wns = 1e7./scan_wvs;
meas_scan_wns = scan_wns;
meas_scan = scan;
clear scan_wvs scan_wns scan;

%Load line positions and amplitudes from lorentzian fit
load(fullfile("..","Saves","HFampfitting.mat"));
amp_fitted_widths(1:4,13:16) = 0.02;
line_positions = amp_fitted_positions;

%Create the axis around the line positions
%coarse axis
coarse = 1e7./linspace(1139, 1141.2, 300);
%fine axis
fine = 1e7./linspace(1139, 1141.2, 20000);
fine(min(abs(fine - line_positions(:))) > 0.05) = [];
scan_wns = sort([coarse fine]);

%Lorentzians
%Positions, Widths, Amps
Lfitting = zeros(size(scan_wns));
bgs = unique(amp_fitted_bgs);
for bg = bgs'
    block_inds = find(amp_fitted_bgs == bg);
    block_positions = line_positions(block_inds);
    block_widths = amp_fitted_widths(block_inds);
    block_amps = amp_fitted_amps(block_inds);

    %For efficiency, focus only on nearby points
    Lparams = [block_positions; block_widths; block_amps; bg];
    mask = min(abs(scan_wns - block_positions)) < 0.05;
    Lfitting(mask) = Lfitting(mask) + multi_lorentz_fun(Lparams, scan_wns(mask), numel(block_amps));
end

%PseudoQuad
%For these, we need to establish the background separate from the signal,
%which means we either fix the background using the above backgrounds or
%inlcude 11 free background parameters in the fit. I will go with the
%former

%NEW way:
%Create a kernel of lorenztians of amplitude 1: a num_data_points x
%num_transitions array. Call it L_kernel. Then the entire spectrum (because
%the strengths are linear, is L_kernel * transition_strengths

%Find the axis to do the fitting on:
m1 = min(abs(meas_scan_wns' - line_positions(:))) < 0.025;
[kernel_axis, sortind] = sort(meas_scan_wns(m1));
meas_scan_trunc = meas_scan(m1);
meas_scan_trunc = meas_scan_trunc(sortind);

% %Temporary - mask out 1/2 to 1/2 transition
% m2 = abs(kernel_axis - 8771.75) >  0.2;
% kernel_axis = kernel_axis(m2);
% meas_scan_trunc = meas_scan_trunc(m2);

%Create the kernel
L_kernel = zeros(numel(kernel_axis), numel(line_positions));
Ext_kernel = zeros(numel(scan_wns), numel(line_positions));
for i = 1:size(L_kernel,2)
    L_kernel(:,i) = multi_lorentz_fun([line_positions(i) amp_fitted_widths(i) 1/amp_fitted_widths(i)], kernel_axis, 1);
    Ext_kernel(:,i) = multi_lorentz_fun([line_positions(i) amp_fitted_widths(i) 1/amp_fitted_widths(i)], scan_wns, 1);
end

%Create the background for the kernel
kernel_bg = zeros(size(kernel_axis));
for i = 1:numel(kernel_axis)
    %use the background of the closest line
    [~,bg_ind] = min(abs(line_positions(:) - kernel_axis(i)));
    kernel_bg(i) = amp_fitted_bgs(bg_ind(1));
end

kb = 0.69503;

%% PseudoQuad Transition Fitting (M1 and PQ)
load(fullfile("..","Saves","PQTransition.mat"));
M1 = zeros(12, 16);
PQ = zeros(12, 16);
for pol = 1:numel(mag_transitions)
    M1 = M1 + abs(mag_transitions{pol}).^2;
    PQ = PQ + abs(elec_transitions{pol}).^2;
end
M1 = M1(:);
PQ = PQ(:);

tempscale = @(T) exp(pretempscale(:)/kb/T );
PQpredicted_spectrum = @(p)  kernel_bg + L_kernel*((p(1)*M1 + p(2)*PQ).*tempscale(p(3)));
PQExt_spectrum = @(p) Ext_kernel*((p(1)*M1 + p(2)*PQ).*tempscale(p(3)));
obj_fun = @(p) meas_scan_trunc - PQpredicted_spectrum(p);

paramsGuess = [0.02 0.09 8];
paramsPQ = lsqnonlin(obj_fun, paramsGuess);

PQ_rel_strength = 100*sum(paramsPQ(2)*PQ)/sum(paramsPQ(1)*M1 + paramsPQ(2)*PQ);

elec_weight = sqrt(paramsPQ(2));
temp = tempscale(paramsPQ(end));
mag_weight = sqrt(paramsPQ(1));
save(fullfile("..","Saves","fittedPQ.mat"), "full_elec_transitions", "elec_weight", "full_mag_transitions", "mag_weight", "temp");

%% Forced Dipole Fitting (M1J and E1J)
load(fullfile("..","Saves","FDTransition.mat"));
M1J = zeros(12,16);
num_elec = size(elec_transitions,1);
for pol = 1:numel(mag_transitions)
    M1J = M1J + abs(mag_transitions{pol}).^2;
end
M1J = M1J(:);

FDpredicted_spectrum = @(p) kernel_bg + L_kernel*((p(1)*M1J...
    + ElecTransition(p(1 + (1:num_elec)), elec_transitions)).*tempscale(p(end)));
FDExt_spectrum = @(p) Ext_kernel*((p(1)*M1J...
    + ElecTransition(p(1 + (1:num_elec)), elec_transitions)).*tempscale(p(end)));
obj_fun = @(p) meas_scan_trunc - FDpredicted_spectrum(p);

paramsGuess = [0.05 0.02*rand(1,num_elec) 8];
paramsFD = lsqnonlin(obj_fun, paramsGuess);

FD_elec_transition = ElecTransition(paramsFD(1 + (1:num_elec)), elec_transitions);
FD_rel_strength = 100*sum(FD_elec_transition)/sum(paramsFD(1)*M1J + FD_elec_transition);

elec_weight = paramsFD(1+(1:num_elec));
temp = tempscale(paramsFD(end));
mag_weight = sqrt(paramsFD(1));
save(fullfile("..","Saves","fittedFD.mat"), "full_elec_transitions", "elec_weight", "full_mag_transitions", "mag_weight", "temp");




%% Plotting
figure("Name","Amplitude Fitting","Units", "inches","Position", [0.5 0.5 3.4 3.3]);

main_ax = axes();
plot(main_ax, 1e7./kernel_axis, meas_scan_trunc, "black");
hold on
% plot(main_ax, 1e7./scan_wns, Lfitting);
Monly = paramsPQ;
Monly(2:end-1) = 0;
plot(main_ax, 1e7./kernel_axis, 0.9*PQpredicted_spectrum(Monly) + 0.1*kernel_bg - 15, "Color", "#db4939");
PQonly = paramsPQ;
PQonly(1) = 0;
plot(main_ax, 1e7./kernel_axis, 4*PQpredicted_spectrum(PQonly) - 3*kernel_bg - 30, "Color", "#3c70c9");
FDonly = paramsFD;
FDonly(1) = 0;
plot(main_ax, 1e7./kernel_axis, 4*FDpredicted_spectrum(FDonly) - 3*kernel_bg - 45, "Color", "#36c781");
% plot(main_ax, 1e7./scan_wns, PQExt_spectrum(paramsPQ));
% plot(main_ax, 1e7./scan_wns, FDExt_spectrum(paramsFD));

legend(["Data" "M1" "PQ" "FD"]);

transitionString = "";
for i = 1:num_elec
    transitionString = transitionString + sprintf("-Path %.1f(%1.0f): %.2f\n", pathJ(i), direction(i), paramsFD(1+i));
end
PQstring = sprintf("PsuedoQuad:\n-Electric Strength: %% %.2f\n-Temp: %.2f", PQ_rel_strength, paramsPQ(end));
FDstring = sprintf("Forced Dipole:\n-Electric Strength %% %.2f\n-Temp %.2f\n" + transitionString, FD_rel_strength, paramsFD(end));
% annotation("textbox", [.2 .5 .3 .3], 'String', PQstring, 'FitBoxToText','on');
% annotation("textbox", [.2 .2 .3 .3], 'String', FDstring, 'FitBoxToText','on');
disp(PQstring)
disp(FDstring)

line_labels = ["1/2 " "3/2 " "5/2 " "7/2 "] + ["1/2" ; "3/2" ; "5/2"];
line_choice = [1 1];
hf_choice = line_positions((4*line_choice(1)-3):(4*line_choice(1)), (4*line_choice(2)-3):(4*line_choice(2)));
%Style xlines by the ground state label (2nd dim)
[~,~,style_labels] = unique(hf_choice(1,:));
[uniquelines, uniqueinds] = unique(hf_choice);
style_labels = repmat(style_labels', 4,1);
styles = ["-" "--" ":"];
for i = uniqueinds'
    xline(1e7/hf_choice(i), "LineStyle",styles(style_labels(i)),"HandleVisibility","off");
end

plot_buffer = 0.02;
xlim(main_ax, [1e7/(max(hf_choice,[],"all")+plot_buffer) 1e7/(min(hf_choice,[],"all")-plot_buffer)]);
xlim(main_ax, [1140.01 1140.03])
ylim(main_ax, [-50 130])
main_ax.Box = "on";
x=0.05;
main_ax.Position = [0.11 (0.59-x) 0.85 (0.375+x)];

main_ax.YLabel.String = sprintf("Fluorescence Amplitude (counts/ms)");
main_ax.XTickLabel([1 end]) = [];
main_ax.XTick([1 end]) = [];

% top_ax = axes();
% top_ax.Position = main_ax.Position;
% top_ax.Color = "none";
% top_ax.XAxisLocation = "top";
% top_ax.YAxisLocation = "right";
% top_ax.XLim = [1e7/main_ax.XLim(2) 1e7/main_ax.XLim(1)];
% top_ax.XDir = "reverse";
% top_ax.YTick = main_ax.YTick;
% top_ax.YTickLabel = [];
% top_ax.YLim = main_ax.YLim;
% top_ax.XLabel.String = "Energy (cm^{-1})";



main_ax.XRuler.TickLabelGapOffset = 0;
main_ax.YRuler.TickLabelGapOffset = 0;
% top_ax.XRuler.TickLabelGapOffset = 0;


%Second plot
main_ax2 = axes();
plot(main_ax2, 1e7./kernel_axis, meas_scan_trunc, "black");
hold on
plot(main_ax2, 1e7./kernel_axis, 0.9*PQpredicted_spectrum(Monly) + 0.1*kernel_bg - 10, "Color", "#db4939");
plot(main_ax2, 1e7./kernel_axis, 8*PQpredicted_spectrum(PQonly) - 7*kernel_bg - 20, "Color", "#3c70c9");
plot(main_ax2, 1e7./kernel_axis, 1.5*FDpredicted_spectrum(FDonly) - 0.5*kernel_bg - 30, "Color", "#36c781");


main_ax2.Box = "on";
main_ax2.Position = [0.11 0.12 0.85 (0.375-x)];
main_ax2.XLim = [1140.197 1140.215];
main_ax2.XLabel.String = "Wavelength (nm)";
main_ax2.XTickLabel(end) = [];
main_ax2.XTick(end) = [];

% top_ax2 = axes();
% top_ax2.Position = main_ax2.Position;
% top_ax2.Color = "none";
% top_ax2.XAxisLocation = "top";
% top_ax2.YAxisLocation = "right";
% top_ax2.XLim = [1e7/main_ax2.XLim(2) 1e7/main_ax2.XLim(1)];
% top_ax2.XDir = "reverse";
% top_ax2.YTick = main_ax2.YTick;
% top_ax2.YTickLabel = [];
% top_ax2.YLim = main_ax2.YLim;

main_ax2.XRuler.TickLabelGapOffset = 0;
main_ax2.YRuler.TickLabelGapOffset = 0;
main_ax2.XTickLabelRotation = 0;
% top_ax2.XRuler.TickLabelGapOffset = 0;

line_choice = [1 2];
hf_choice = line_positions((4*line_choice(1)-3):(4*line_choice(1)), (4*line_choice(2)-3):(4*line_choice(2)));
%Style xlines by the ground state label (2nd dim)
[~,~,style_labels] = unique(hf_choice(1,:));
[uniquelines, uniqueinds] = unique(hf_choice);
style_labels = repmat(style_labels', 4,1);
styles = ["-" "--" ":"];
for i = uniqueinds'
    xline(main_ax2, 1e7/hf_choice(i), "LineStyle",styles(style_labels(i)),"HandleVisibility","off");
end

set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','FontName'),'FontName',"Times New Roman");







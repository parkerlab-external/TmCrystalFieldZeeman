load("MagnetResults.mat");
% B_measured = [155, 130, 100, 98, 85];

sample_B = linspace(70, 160, 20);
sample_B = sample_B/10000*1000;
% sample_B = 0;
linewidths = linspace(0.004, 0.005, 3);
linewidths = 0.0045;

figure("Name","Polarized Spectra","Units", "inches","Position", [0.5 0.5 3.4 3.3]);
main_ax = axes();



xaxisinput = 1e7./[probe1 probe2] - 80e6/29979245800;

hold on;

manual_bg= -0.6;
num_rand_samples = 40;
theory_probe_ratios = zeros(numel(linewidths), numel(sample_B));
theory_probe_ratio_stds = zeros(numel(linewidths), numel(sample_B));
for i = 1:numel(linewidths)
    linewidth = linewidths(i);
    for j = 1:numel(sample_B)
        probe_samples = zeros(2,num_rand_samples);
        probe_ratio_samples = zeros(1,num_rand_samples);
        for k = 1:num_rand_samples
            [xaxis, temp] = FullInteraction(sample_B(j)/1000, linewidth, "unpolarized1","xaxis", xaxisinput, "N", 50);
%             probe_samples(1,k) = interp1(1e7./xaxis, temp, probe1);
            probe_samples(1,k) = temp(1)+manual_bg;
%             probe_samples(2,k) = interp1(1e7./xaxis, temp, probe2);
            probe_samples(2,k) = temp(2)+manual_bg;
            probe_ratio_samples(k) = probe_samples(1,k)/probe_samples(2,k);
        end
        A = mean(probe_samples(1,:));
        dA = std(probe_samples(1,:))/sqrt(num_rand_samples);
        B = mean(probe_samples(2,:));
        dB = std(probe_samples(2,:))/sqrt(num_rand_samples);
        theory_probe_ratios(i,j) = A/B;
        theory_probe_ratio_stds(i,j) = A/B * sqrt((dA/A)^2 + (dB/B)^2);
        theory_probe_ratios(i,j) = mean(probe_ratio_samples);
        theory_probe_ratio_stds(i,j) = std(probe_ratio_samples);
        disp(j)
    end
end

for i = 1:size(theory_probe_ratios,1)
%     errorbar(sample_B, theory_probe_ratios(i,:), theory_probe_ratio_stds(i,:), 'o');
    temp = plot(sample_B, theory_probe_ratios(i,:), "Color", "#db4939");
    hold on
    fillx = [sample_B sample_B(end:-1:1)];
    filly = [(theory_probe_ratios(i,:)-theory_probe_ratio_stds(i,:)) (theory_probe_ratios(i,end:-1:1)+theory_probe_ratio_stds(i,end:-1:1))];
    fill(fillx, filly, temp.Color, "FaceAlpha",0.3, 'EdgeColor','none');
end

temp = errorbar(B_measured/10, probe_ratios, probe_ratio_stds, 'ok');
errorbar(B_measured/10, probe_ratios, B_err/10, "horizontal", "Color", temp.Color, "LineStyle","none");

main_ax.XLabel.String = "Magnetic Field (mT)";
main_ax.YLabel.String = "Fluorescence Ratio";
main_ax.XLim = [7 16];


set(findall(gcf, "-property", "FontSize"), "FontSize", 10);
set(findall(gcf, "-property", "FontName"), "FontName", "Times New Roman");
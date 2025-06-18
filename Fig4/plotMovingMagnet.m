%Estimating magnetic field
dists = [93 98 105 109 113]; %mm, between black pieces that hold the magnets
%B_measured = [147.44 120.32 101.32 90.54 81.43];

%B_measured is measured in Gauss, second dimension is 0.1in
B_measured = [194.02 172.51 157.78 153.21 147.44 149.53 162.39 187.23 225.15;
    201.45 168.26 146.22 131.90 124.33 120.32 120.37 123.91 131.47;
    143.80 126.23 114.20 107.02 102.41 101.32 101.54 105.42 112.26;
    128.92 113.69 102.86 95.74 91.83 90.54 91.29 94.51 99.12;
    112.09 99.36 91.08 86.03 82.35 81.38 81.43 83.79 88.30];

%Fit to field 3z^2 - 3(x^2 + y^2)/2 + C, measured in z so just fit to
%parabola
figure("Name", "B polyfitting")
z = (1:size(B_measured, 2))*2.54; %mm scale
fit_z = linspace(min(z), max(z), 100);
B_min = zeros(1, size(B_measured,1));
B_err = zeros(1, size(B_measured,1));
for i = 1:numel(B_min)
    B_curve = B_measured(i, :);
    p = polyfit(z, B_curve, 2);

    %a(x-b)^2 + c = ax^2 - 2ab*x + ab^2 + c
    small_b = p(2)/(-2*p(1));
    small_c = p(3) - p(1)*small_b^2;
    
    %3z^2 - 3(x^2 + y^2)/2
    B_func = @(z, r) p(1)*((z-small_b).^2 - 0.5*r.^2);
    B_mag = @(z, r) sqrt(small_c^2 + 2*small_c*B_func(z, r));

    B_min(i) = B_mag(small_b, 0);
    rand_samples = randn([2, 1000])*2;
    B_err(i) = std(B_mag(rand_samples(1,:) + small_b, rand_samples(2,:)));


    plot(z, B_curve)
    hold on
    plot(fit_z, polyval(p,fit_z));
    plot(fit_z, B_func(fit_z, 0) + small_c);
    plot(fit_z, B_mag(fit_z, 0));
    errorbar(small_b, B_min(i), B_err(i));
end
hold off

B_measured = B_min;






%Loading data
files = dir("ForFigure\*.mat");

figure("Name","MovingMagnet");

main_ax = axes();
hold on

for i = 1:numel(files)
    file = files(i);
    load(fullfile(file.folder,file.name));
    mean_data = squeeze(mean(data_out.raw_data,2));
    data_error = squeeze(std(data_out.raw_data,[],2))/sqrt(size(data_out.raw_data,2));
    decay_data = mean_data(:,22:end);
    decay_errors = data_error(:,22:end);
%     [scan1, tc1, ~, paramsVar1] = decay_fit(1:size(decay_data,2), decay_data);
    [scan, tc, ~, paramsVar] = decay_fit(1:size(decay_data,2), decay_data, "decay_errors",decay_errors);
    wv_samples = data_out.wv_samples;
    scan_std = paramsVar(1:numel(scan));
    
    scan_wvs = interp1(find(wv_samples ~= 0),wv_samples(wv_samples ~= 0), 1:numel(wv_samples));
    errorbar(main_ax, scan_wvs, scan, paramsVar(1:numel(scan)));

    %Put everything on the same axis
    if i == 1
        fixed_wvs = scan_wvs;
        fixed_scans = zeros(numel(files), numel(scan));
        fixed_stds = zeros(numel(files), numel(scan));
        fixed_scans(i,:) = scan;
        fixed_stds(i,:) = scan_std;
    else
        fixed_scans(i,:) = interp1(scan_wvs, scan, fixed_wvs, 'linear', 'extrap');
        fixed_stds(i,:) = interp1(scan_wvs, scan_std, fixed_wvs, 'linear', 'extrap');
    end
end

%Find a ratiometric measurement on each iso-B
ratios = zeros(numel(files), numel(fixed_wvs), numel(fixed_wvs));
ratio_stds = zeros(size(ratios));
for i = 1:size(ratios,1)
    A = fixed_scans(i,:);
    B = fixed_scans(i,:)';
    dA = fixed_stds(i,:);
    dB = fixed_stds(i,:)';
    ratios(i,:,:) = A./B;
    ratio_stds(i,:,:) = (A./B).*sqrt((dA./A).^2 + (dB./B).^2);
end

%from this we find the two wavelengths for the comparison measurement
temp = squeeze(ratios(1,:,:) - ratios(5,:,:))./squeeze(sqrt(ratio_stds(1,:,:).^2 + ratio_stds(5,:,:).^2));
[probei, probej] = find(temp == max(temp, [], "all"));
probei = 134;
probej = 154;
probe1 = fixed_wvs(probei);
probe2 = fixed_wvs(probej);

figure;
plot(fixed_wvs, fixed_scans);
xline(probe1);
xline(probe2,"red");

%Finally confine the measurements to the probes
probe_ratios = zeros(1,numel(files));
probe_ratio_stds = zeros(1,numel(files));
for i = 1:numel(probe_ratios)
    probe_ratios(i) = ratios(i,probej, probei);
    probe_ratio_stds(i) = ratio_stds(i,probej, probei);
end

figure("Name", "Ratio Measurement vs B");
errorbar(B_measured, probe_ratios, probe_ratio_stds);

save("MagnetResults.mat", "B_measured", "B_err", "probe_ratios", "probe_ratio_stds", "probe1", "probe2");

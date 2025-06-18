%No Field
load(fullfile("Important Magnet Data","vspectrum-range12024-11-12-1655.mat"));
wv_samples = data_out.wv_samples;
interp_wvs = interp1(find(wv_samples ~= 0), wv_samples(wv_samples ~= 0), 1:numel(wv_samples));
mask = data_out.set_sms < 1.2e4 & ([0 ; abs(diff(interp_wvs'))] < 6e-5) & interp_wvs' > 1140.01 & interp_wvs' < 1140.03;
mean_data = squeeze(mean(data_out.raw_data(mask,:,:),2));

%Additional removal of overlap
scan_wvs1 = interp_wvs(mask);
scan_wvs1(362:623) = [];

decay_data = mean_data(:,5:end);
decay_data(362:623, :) = [];


[scan1, tc1, bg1, paramVars1] = decay_fit(1:size(decay_data,2), decay_data);
[scan_wvs1, sortind] = sort(scan_wvs1);
scan1 = scan1(sortind);

figure("Name", "BFieldPolarized");
main_axis = axes();
plot(main_axis, scan_wvs1, scan1);
hold on;

%One polarization
load(fullfile("Important Magnet Data","vspectrum2024-11-12-0850.mat"));
wv_samples = data_out.wv_samples;
interp_wvs = interp1(find(wv_samples ~= 0), wv_samples(wv_samples ~= 0), 1:numel(wv_samples));
mask = data_out.set_sms < 1.2e4 & ([0 ; abs(diff(interp_wvs'))] < 6e-5) & interp_wvs' > 1140.01 & interp_wvs' < 1140.03;
mean_data = squeeze(mean(data_out.raw_data(mask,:,:),2));
decay_data = mean_data(:,5:end);

[scan2, tc2, bg2, paramVars2] = decay_fit(1:size(decay_data,2), decay_data);
scan_wvs2 = interp_wvs(mask);

plot(main_axis, scan_wvs2, scan2)

%Other polarization
load(fullfile("Important Magnet Data","vspectrum2024WP-45-11-11-2058.mat"));
wv_samples = data_out.wv_samples;
interp_wvs = interp1(find(wv_samples ~= 0), wv_samples(wv_samples ~= 0), 1:numel(wv_samples));
mask = data_out.set_sms < 1.2e4 & ([0 ; abs(diff(interp_wvs'))] < 6e-5) & interp_wvs' > 1140.01 & interp_wvs' < 1140.03;
mean_data = squeeze(mean(data_out.raw_data(mask,:,:),2));
decay_data = mean_data(:,5:end);

[scan3, tc3, bg3, paramVars3] = decay_fit(1:size(decay_data,2), decay_data);
scan_wvs3 = interp_wvs(mask);

plot(main_axis, scan_wvs3, scan3)

save("PolarizedSpectra", "scan_wvs1", "scan1", "scan_wvs2", "scan2", "scan_wvs3", "scan3");
load('vspectrum2024-10-08-1717.mat'); %Magnet
% load('vspectrum2024-10-08-1553.mat'); %No magnet
mean_data = squeeze(mean(data_out.raw_data,2));
decay_data = mean_data(:,17:end);
scan = decay_fit(1:size(decay_data,2), decay_data);
wv_samples = data_out.wv_samples;
scan_wvs = interp1(find(wv_samples ~=0), wv_samples(wv_samples~=0), 1:numel(wv_samples));

plot(scan_wvs, scan);


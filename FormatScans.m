%% Process raw spectrum of choice to be used for the position fitting, amplitude fitting
% This script will change depending on what data set you want to use, but
% should always output into the Spectra_formatted folder in the required
% format

%% Site I
%load preprocessed scan
load(fullfile("InitialSpectra","SiteISpectrum.mat"))
scan_wvs = reshape(scan_wvs, [], 1);
scan = reshape(scan, [], 1);
scan_err = reshape(scan_err, [], 1);

%Already has been shifted by 80 MHz

%% Site II
% load(fullfile("InitialSpectra", "vspectrum2024-10-02-2054.mat"));
% N = numel(data_out.wv_samples);
% for i = 1:N
%     set_sms(i) = data_out.set_operations{i}(end).meas_array(4,end);
% end
% 
% mask = data_out.set_results' == "success" & set_sms < 2.3e4;
% 
% %Statistics
% mean_data = squeeze(mean(data_out.raw_data(mask,:,:),2));
% std_data = squeeze(std(data_out.raw_data(mask,:,:), 0, 2));
% 
% %Preallocate
% num_points = size(mean_data, 1);
% pump = zeros(1,num_points);
% probe = zeros(1,num_points);
% pump_probe = zeros(1, num_points);
% 
% decay_axis = 0:81;
% for i = 1:size(mean_data,1)
%     pump_decay = mean_data(i,17:98);
%     pump_probe_decay = mean_data(i,125:206);
%     probe_decay = mean_data(i,229:310);
% 
%     [pump(i),tc,bg] = decay_fit(decay_axis, pump_decay);
%     [pump_probe(i),~,~,pErr] = decay_fit(decay_axis, pump_probe_decay, 'FixedBackground',bg, "tc_lb", 13, "tc_ub", 30,"decay_errors", std_data(i,125:206));
%     pump_probe_err(i) = pErr(1);
%     [probe(i),~,~,pErr] = decay_fit(decay_axis, probe_decay, 'FixedBackground',bg, "tc_lb", 13, "tc_ub", 30, "decay_errors", std_data(i, 229:310));
%     probe_err(i) = pErr(1);
% end
% 
% scan_wvs = data_out.wv_samples(mask);
% scan = pump + probe - pump_probe;
% scan_wvs = reshape(scan_wvs, [], 1);
% scan = reshape(scan, [], 1);
% scan_err = sqrt(probe_err.^2 + pump_probe_err.^2);
% scan_err = reshape(scan_err, [], 1);
% 
% c = 299792458;
% %Measuring the +1 order of the AOM, shifted by 80 MHz, so we need to
% %subtract 80MHz to get the real wavelength
% scan_wvs =  1./(1./scan_wvs - (80e6)/(c*(1e9)) );


%% Save

save(fullfile("Saves","rough_position_scan.mat"), "scan_wvs","scan", "scan_err");

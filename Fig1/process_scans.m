pumpprobe_files = [
    ".\vspectrum2024-09-30-1735.mat";
    ".\vspectrum2024-10-02-2054.mat";
    ];

for j = 1:numel(pumpprobe_files)
    filename = pumpprobe_files(j);
    load(filename);

    mean_data = squeeze(mean(data_out.raw_data,2));
    set_sms = zeros(size(data_out.wv_samples));
    for i = 1:numel(set_sms)
        lastop = data_out.set_operations{i}(end);
        set_sms(i) = lastop.meas_array(lastop.sms_ind,end);
    end
    mask = set_sms < 1.5e4 & data_out.set_results == "success";
    mean_data = mean_data(mask,:);
    
    pump_decay = mean_data(:,16:102);
    pumpprobe_decay = mean_data(:,124:207);
    probe_decay = mean_data(:,228:end);
    
    pump_tcs = zeros(1, size(pump_decay, 1));
    pumpprobe_tcs = zeros(1, size(pumpprobe_decay, 1));
    probe_tcs = zeros(1, size(probe_decay,1));
    
    pump = zeros(size(mean_data,1),1);
    probe = zeros(size(pump));
    pumpprobe = zeros(size(pump));
    
    bg = 10;
    
    for a = 1:size(mean_data,1)
        [pump(a), pump_tcs(a)] = decay_fit(0:size(pump_decay,2)-1, pump_decay(a,:), FixedBackground=bg);
        [pumpprobe(a), pumpprobe_tcs(a)] = decay_fit(0:size(pumpprobe_decay,2)-1, pumpprobe_decay(a,:), FixedBackground=bg);
        [probe(a), probe_tcs(a)] = decay_fit(0:size(probe_decay,2)-1, probe_decay(a,:), FixedBackground=bg);
    end
    
    pump_wv = shift_wavelengths(mean(data_out.interp_pump_wavelengths));
    scan_wvs = shift_wavelengths(data_out.wv_samples(mask));

    save(append(erase(filename, ".mat"),"_processed.mat"), "data_out", "pump", "pump_tcs", "probe", "probe_tcs", "pumpprobe", "pumpprobe_tcs", "pump_wv", "scan_wvs");
end
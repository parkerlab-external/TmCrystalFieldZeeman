files = dir("ForFigure\*.mat");

%get the background and time constant:
% load("ForFigure\vspectrum2024-10-11-1344.mat");
% mean_data = squeeze(mean(data_out.raw_data(:,:,23:end),2))/2;
% [~,tc,bg] = decay_fit(1:size(mean_data,2), mean_data);


load("MagnetResults.mat");

%empty arrays
R = zeros(10,numel(files));
p1 = R;
p2 = R;
Rerror_scaling = zeros(10, numel(files));
p1error_scaling = Rerror_scaling;
p2error_scaling = Rerror_scaling;



num_repeats = 10;
for i = 1:numel(files)
    load(fullfile("ForFigure/",files(i).name));
    %Find the correct probe wavelength
    m = data_out.wv_samples ~= 0;
    scan_wvs = interp1(find(m), data_out.wv_samples(m), 1:numel(data_out.wv_samples));
    [~,probei] = min(abs(scan_wvs - probe1));
    [~,probej] = min(abs(scan_wvs - probe2));

    for j = 3:num_repeats
        raw_data = data_out.raw_data(:,1:j,23:end);
        mean_data = squeeze(mean(raw_data,2));
        std_data = squeeze(std(raw_data,[],2))/sqrt(j);

        %probe decay
        [amps, ~, ~, paramVar] = decay_fit(1:size(mean_data,2), mean_data, "decay_errors",std_data);
        probe1err = sqrt(interp1(scan_wvs, paramVar(1:numel(amps)).^2, probe1));
        probe1amp = interp1(scan_wvs, amps, probe1);
        probe2err = sqrt(interp1(scan_wvs, paramVar(1:numel(amps)).^2, probe2));
        probe2amp = interp1(scan_wvs, amps, probe2);

%         figure(1);
%         plot(scan_wvs, amps);
%         hold on

        %Save in the empties
        p1(j,i) = probe1amp;
        p2(j,i) = probe2amp;
        R(j,i) = probe1amp/probe2amp;
        p1error_scaling(j,i) = probe1err;
        p2error_scaling(j,i) = probe2err;
        Rerror_scaling(j,i) = R(j,i)*sqrt((probe1err/probe1amp)^2 + (probe2err/probe2amp)^2);
    end
end

%% Find the average slope to get an idea of the R per G sensitivity
p = polyfit(B_measured/10, R(10,:), 1);
best_sens = Rerror_scaling/p(1);
best_sens = best_sens(3:10,:);
integration_time = (3:10)*0.2;

dimless_sens = zeros(1,5);
for i = 1:numel(dimless_sens)
    ptemp = polyfit(1./sqrt(2*integration_time), best_sens(:,i), 1);
    dimless_sens(i) = ptemp(1);
end




%% Extra
% %How does the error scale with signal?
% figure;
% for i = 1:numel(files)
%     load(fullfile("ForFigure/", files(i).name));
%     raw_data = data_out.raw_data(:,:,23:end)/2;
% 
%     mean_data = squeeze(mean(raw_data, 2));
%     std_data = squeeze(std(raw_data,0,2));
% 
%     scatter(log(mean_data(:)), log(std_data(:)), '.');
%     hold on;
% end
% xlabel("log(Signal)")
% ylabel("log(Noise)");
% 
% x = linspace(3, 6);
% plot(x, 0.5*x);
% 
% %How does the error scale with integration time?
% inverseRootHertz = [];
% signalSTD = [];
% figure;
% for i = 1:numel(files)
%     load(fullfile("ForFigure/", files(i).name));
%     raw_data = data_out.raw_data(:,:,23:end);
%     for j = 2:size(raw_data,2)
%         mean_data = squeeze(mean(raw_data(:,1:j,:),2));
%         std_data = squeeze(std(raw_data(:,1:j,:),0,2))/sqrt(j-1);
% 
%         [amps, ~, ~, paramsVar] = decay_fit(1:size(mean_data,2), mean_data, "decay_errors", std_data, "FixedBackground", bg, "tc_lb", tc, "tc_ub", tc);
%         scatter(1/sqrt(j*0.1), mean(paramsVar(1:numel(amps))), ".k");
%         inverseRootHertz(end+1) = 1/sqrt(j*0.1);
%         signalSTD(end+1) = mean(paramsVar(1:numel(amps)));
%         hold on;
%     end
% end


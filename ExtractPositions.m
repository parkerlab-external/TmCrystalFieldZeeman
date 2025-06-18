function ExtractPositions(estimates, widths, num_curves)
    %Provided with estimates of lines and general widths, this function
    %will generate the true line centers, provided the estimates are decent

    %estimates: 1xN array of line position estimates, where N must be the
    %number of lines (this function does not combine two similar lines into
    %one)

    %widths: 1xN array of line widths, the range in which the lines will be
    %fit with


    %Establish the regions where the fitting will be done
    regions = [(estimates - widths/2) ; (estimates + widths/2)]; % 2xN array of ranges, may be overlapped
    regions = shared_regions(regions); %Remove overlaps to prevent double counting

    %Load the processed scan
    load(fullfile("Saves","rough_position_scan.mat"));
    scan_wns = 1e7./scan_wvs;

    %Prepare amplitude and position arrays for the lines
    line_amps = zeros(size(estimates));
    line_pos = zeros(size(estimates));
    line_widths = zeros(size(estimates));
    line_pos_err = line_pos;

    for region_num = 1:size(regions,2)
        region = regions(:,region_num);
        region = sort(region);
        m = scan_wns > region(1) & scan_wns < region(2);
        focus_scan = scan(m);
        focus_err = scan_err(m);
        focus_wns = scan_wns(m);

        if ~isempty(focus_scan)
            [params, ~, paramErr] = lorentzian_fit(focus_wns, focus_scan, num_curves, "NumRandSamples", 20, "scan_err", focus_err);

            %Relevant parameters
            positions = params(1:num_curves);
            positions_err = paramErr(1:num_curves);
            widths = params(num_curves+1:2*num_curves);
            amps = params(2*num_curves+1:3*num_curves);
            amps = amps.*widths;

            figure;
            scatter(focus_wns, focus_scan)
            hold on
            plot(focus_wns, multi_lorentz_fun(params, focus_wns, num_curves));
            xline(positions)
            errorbar(positions, mean(focus_scan), positions_err, "horizontal", "LineStyle","none")

            %Modify parameters for the closest estimated line:
            for j = 1:num_curves
                pos = positions(j);
                pos_err = positions_err(j);
                [~, min_ind] = min(abs(estimates - pos));
                
                line_pos(min_ind) = (line_amps(min_ind)*line_pos(min_ind) + amps(j)*pos)/(line_amps(min_ind) + amps(j)); %Modify weighted average position
                line_pos_err(min_ind) = sqrt( (line_amps(min_ind)*line_pos_err(min_ind))^2 + (amps(j)*pos_err)^2)/(line_amps(min_ind) + amps(j));
                line_amps(min_ind) = line_amps(min_ind) + amps(j);
                if num_curves == 1
                    line_widths(min_ind) = widths(j);
                end
            end
            xline(line_pos(min_ind), "blue")
            errorbar(line_pos(min_ind), 2/3*(max(focus_scan) - min(focus_scan)) + min(focus_scan), line_pos_err(min_ind), "horizontal", "LineStyle","none");
        end
    end

    line_positions = line_pos;
    line_positions_err = line_pos_err;
    line_amplitudes = line_amps;
    save(fullfile("Saves","coarse_lines.mat"), "line_positions", "line_amplitudes", "line_positions_err", "line_widths");
end
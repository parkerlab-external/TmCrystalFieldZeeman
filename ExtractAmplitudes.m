%This script takes positions as a given, and forces lorentzian fits on the
%given positions. This is designed to give amplitudes for hyperfine levels

load(fullfile("Saves","rough_position_scan.mat"));
scan_wns = 1e7./scan_wvs;

%Don't use the crystal field coarse levels, use the ones from the level
%fitting (or from the original COM fitting, whichever is better):
load(fullfile("Saves","coarse_levels.mat"))
coarse_positions = repelem(sort(eE),4)' - repelem(sort(gE),4);
load(fullfile("Saves", "CFHFLevels.mat"));
HF_shifts = HFShifts(17:28) - HFShifts(1:16)';
mFsource = repmat(mF(17:28), [1 16]);
mFdest = repmat(mF(1:16)', [12 1]);


%Set the "look range" to be the maximum hyperfine splitting
width = 0.4*max(abs(HF_shifts),[],"all");

line_positions = coarse_positions + HF_shifts;
regions = [line_positions(:)- width, line_positions(:) + width];
regions = shared_regions(regions'); %Remove overlaps

%Preallocating
amp_fitted_shifts = zeros(size(line_positions));
amp_fitted_widths = ones(size(line_positions));
amp_fitted_amps = zeros(size(line_positions));
amp_fitted_bgs = zeros(size(line_positions));
amp_fitted_position_errs = zeros(size(line_positions));
amp_fitted_positions = zeros(size(line_positions));

for region_num = 1:size(regions,2)
    region = regions(:,region_num);
    region = sort(region);

    m = scan_wns > region(1) & scan_wns < region(2);
    linemask = find(line_positions > region(1) & line_positions < region(2));

    focus_scan = scan(m);
    focus_err = scan_err(m);
    focus_wns = scan_wns(m);
    
    if ~isempty(focus_scan)
        figure;

        %Extract unique line positions to know how many curves to fit
        [focus_lines, focus_inds] = uniquetol(line_positions(linemask),1e-10);
        line_tolerance = 0.7*mean(abs(diff(sort(focus_lines))));

        %The objective function takes an offset, width, amplitude array,
        %and background as parameters
        obj_fun = @(params) focus_scan - fixed_lorentz_fun(focus_wns, focus_lines, params);

        %Construct bounds and guesses for params: [amps x num_fixed offset width bg]
        lb = [zeros(1,numel(focus_lines)) -1*mean(abs(diff(sort(focus_lines))))/3  mean(diff(sort(focus_wns))) 0];
        ub = [1.1*max(focus_scan)*ones(1,numel(focus_lines)) 1*mean(abs(diff(focus_lines)))/3  (max(focus_wns)-min(focus_wns)) 1.1*max(focus_scan)];
        paramGuess = lb + rand(size(lb)).*(ub - lb);

        [params, ~, ~, ~, ~, ~, J] = lsqnonlin(obj_fun, paramGuess, lb, ub);

        T = inv(J'*J);
        paramsErr = abs(T*J.'*diag(focus_err(:)));

        %Block (linemask) properties
        amp_fitted_bgs(linemask) = params(end);
        amp_fitted_widths(linemask) = params(end-1);
        amp_fitted_shifts(linemask) = params(end-2);
        amp_fitted_position_errs(linemask) = paramsErr(end-2);

        %Line Specific properties
        amp_fitted_positions(linemask) = line_positions(linemask) + params(end-2);
        amp_fitted_amps(linemask(focus_inds)) = params(1:numel(focus_lines))*params(end-1);

        %Plot
        nlines = numel(focus_lines);
        plot_params = [amp_fitted_positions(linemask(focus_inds))' repelem(params(end-1), nlines) amp_fitted_amps(linemask(focus_inds))'/params(end-1) params(end)];

        plot(1e7./focus_wns, focus_scan)
        hold on
        plot(1e7./focus_wns, multi_lorentz_fun(plot_params, focus_wns, nlines))
        xline(1e7./focus_lines, "black")
        xline(1e7./(focus_lines + params(end-2)), "blue");
        text(min(xlim), max(ylim), sprintf("FWHM: %.4f", 2*mean(params(end-1))), 'Horiz','left', 'Vert','top');
    end
end

save(fullfile("Saves","HFampfitting.mat"), "amp_fitted_shifts", "amp_fitted_widths","amp_fitted_amps","amp_fitted_bgs", "amp_fitted_positions", ...
  "mFsource", "mFdest", "amp_fitted_position_errs");



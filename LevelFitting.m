function LevelFitting(line_positions, line_positions_err)
measured_lines = line_positions./line_positions_err;

%% Fitting
NUM_MISSING = 1;
NUM_GROUND = 4;
NUM_EXCITED = 3;
G_labels = compose("%c", 'A':char('A' + NUM_GROUND-1)); %["A" "B" "C" "D"]
E_labels = compose("%c", '1':char('1' + NUM_EXCITED-1)); %["1" "2" "3"]
NUM_MEASURED = numel(measured_lines);

[temp1, temp2] = ndgrid(G_labels, E_labels);
label_strings = strcat(temp1, temp2)';
clear temp1 temp2

order_cands = 1:numel(label_strings);

full_orders = zeros(size(order_cands));
for gState = 1:NUM_GROUND
    %Assign the first missing element
    for i = 1:size(order_cands,1)
        full_orders(i,order_cands(i,1)) = (gState-1)*NUM_EXCITED + 1;
    end
    order_cands(:,1) = [];

    %Remaining labels
    states_remaining = (gState-1)*NUM_EXCITED + (2:NUM_EXCITED);


    orig_size = size(full_orders,1);
    sample_factor = nchoosek(numel(order_cands(i,1:(end-NUM_GROUND+gState))), numel(states_remaining));

    %Create new arrays
    full_orders = repelem(full_orders, sample_factor, 1);
    new_order_cands = zeros(orig_size*sample_factor, size(order_cands,2)-numel(states_remaining));

    for i = 1:orig_size
        %Sample the remaining candidate indices
        samples = nchoosek(order_cands(i,1:(end-NUM_GROUND+gState)), numel(states_remaining));

        for j = 1:sample_factor
            ind = (i-1)*sample_factor + j;
            full_orders(ind, samples(j,:)) = states_remaining;
            new_order_cands(ind, :) = sort(find(full_orders(ind,:) == 0));
        end
    end
    order_cands = new_order_cands;
end

clear new_order_cands order_cands orig_size sample_factor samples states_remaining gState ind

G_matrix = repelem(-1*eye(NUM_GROUND),NUM_EXCITED,1);
E_matrix = 1*eye(NUM_EXCITED);
E_matrix = repmat(E_matrix, NUM_GROUND, 1);
full_coeff_matrix = [G_matrix E_matrix];
clear G_matrix E_matrix
%Not fitting to the ground state lowest level "A1"
full_coeff_matrix(:,1) = [];

delete_labels = 1:numel(label_strings);
delete_labels = nchoosek(delete_labels, NUM_MISSING);


total_checked = 0;
hold_top = 10;
best_orders = zeros(hold_top, numel(label_strings)-NUM_MISSING);
best_full_orders = zeros(hold_top, numel(label_strings));
best_spacing = zeros(hold_top, size(full_coeff_matrix,2));
best_residuals = Inf*(ones(hold_top,1));
for i = 1:size(delete_labels, 1)
    %delete_labels removes a certain index from the coefficient matrix,
    %meaning it removes a certain LABEL from the ordered states, so we
    %can't just remove that same index from the ordered state because it
    %might not be associated with the correct label to be deleted

    missing_coeff_matrix = full_coeff_matrix;
    missing_coeff_matrix(delete_labels(i,:),:) = [];
    
    %Loop through each ordering and remove the LABELS specified by
    %delete_labels
    for j = 1:size(full_orders, 1)
        missing_orders = full_orders(j,:);
        missing_orders(ismember(missing_orders, delete_labels(i,:))) = [];
        
        %Construct ordered measured values
        missing_measured = zeros(numel(label_strings),1);
        missing_measured(missing_orders) = measured_lines;
        missing_err = missing_measured;
        missing_err(missing_orders) = measured_lines./line_positions;
        if missing_measured(delete_labels(i,:)) ~= 0
            disp("OK")
        end
        missing_measured(delete_labels(i,:)) = [];
        missing_err(delete_labels(i,:)) = [];
        


        level_spacing = (missing_coeff_matrix.*missing_err)\missing_measured;
        residual = norm((missing_coeff_matrix.*missing_err)*level_spacing - missing_measured);

        if residual < max(best_residuals)
            if ~ismember(missing_orders, best_orders, 'rows')
                best_residuals(end) = residual;
                best_orders(end,:) = missing_orders;
                best_full_orders(end,:) = full_orders(j,:);
                best_spacing(end,:) = level_spacing;
    
                [best_residuals, sortind] = sort(best_residuals);
                best_orders = best_orders(sortind,:);
                best_full_orders = best_full_orders(sortind,:);
                best_spacing = best_spacing(sortind,:);
            end
        end
        total_checked = total_checked+1;
    end
end
clear residual missing_orders level_spacing
clear i j missing_measured missing_coeff_matrix

%% Plotting
load(fullfile("Saves", "rough_position_scan.mat"))
plot_ind = 1;
full_predicted = full_coeff_matrix*best_spacing(plot_ind,:)';
figure("Name", "Level Fitting");
plot((1e7)./scan_wvs, scan)
xline(measured_lines.*line_positions_err)
xline(full_predicted, "--r", label_strings(:))

%% Important Saves
%Levels
[full_spacing, sind] = sort([0 best_spacing(1,:)]);
shifted_spacing = full_spacing - mean(full_spacing(1:NUM_GROUND));
%Create basis transform which does the above
shift_T = zeros(numel(full_spacing), numel(full_spacing));
shift_T(:,1:NUM_GROUND) = -1/NUM_GROUND; %Average of ground
shift_T = shift_T + eye(size(shift_T));

gE = shifted_spacing(1:NUM_GROUND);
eE = shifted_spacing((1:NUM_EXCITED) + NUM_GROUND);

A = full_coeff_matrix(best_orders(1,:), :);
T = inv(A.'*A);
% level_Cov = T*A.'*line_positions_err'*line_positions_err*A*T.';
level_Cov = T*A.'*diag(line_positions_err(:).^2)*A*T.';
level_Cov = padarray(level_Cov, [1 1], 0, 'pre');
level_Cov = level_Cov(sind, sind);
level_Cov = shift_T * level_Cov * shift_T';
level_err = sqrt(diag(level_Cov));

gE_err = level_err(1:NUM_GROUND);
eE_err = level_err((1:NUM_EXCITED) + NUM_GROUND);

save(fullfile("Saves","coarse_levels.mat"), "gE", "eE", "gE_err", "eE_err");

%Amplitudes with level-fitted positions
line_positions = eE' - gE;
% loaded_amps = line_amplitudes;
% line_amplitudes = zeros(size(line_positions));
% line_amplitudes(best_orders(1,:)) = loaded_amps;
% line_amplitudes = line_amplitudes(1:end, end:-1:1);

% save("Lines\level_fitted_lines", "line_positions", "line_amplitudes");

%Amplitudes with original positions
line_positions = line_positions(1:end, end:-1:1);
line_positions(best_orders(1,:)) = measured_lines;
line_positions = line_positions(1:end, end:-1:1);
% save("Lines\level_raw_lines", "line_positions", "line_amplitudes");
end





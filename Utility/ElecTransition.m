function total_transition = ElecTransition(params, full_elec_transitions)
    %This function coherently adds transitions and returns the
    %total_transitions

    sum_elec_transitions = cell(1, size(full_elec_transitions,2));
    
    %Coherently add paths
    for pol = 1:size(full_elec_transitions,2)
        sum_elec_transitions{pol} = zeros(size(full_elec_transitions{1}));
        for i = 1:numel(params)
            sum_elec_transitions{pol} = sum_elec_transitions{pol} + params(i)*full_elec_transitions{i, pol};
        end
    end

    %Square elements and add up the polarizations
    total_transition = zeros(size(sum_elec_transitions{1}));
    for pol = 1:numel(sum_elec_transitions)
        total_transition = total_transition + abs(sum_elec_transitions{pol}).^2;
    end

    %Elongate
    total_transition = total_transition(:);

end
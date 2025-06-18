function shared_regions = shared_regions(regions)
%This function takes regions in the 2xN matrix form as
% [ x_0min x_1min ... x_Nmin ]
% [ x_0max x_1max ... x_Nmax ]

%and returns ranges with no overlap
% [x_0min ... x_Mmin ]
% [x_0min ... x_Mmax ]

%The input regions do not have to be sorted in any way
%The output regions are sorted from lowest to highest

%METHOD: assign each pair to a "bit", and sort all points of interest,
%keeping track of which points belong to which pairs. The minimums
%correspond to a positive bit flip, and maximums to a negative. So a range
%is complete when all pairs in the range have a partner, or all bits are 0.
%So we can do a trick by assigned positive and negative integers
%the minimums and maximums of each range and tracking the sums through the
%sorted points

%Isolate minumums and maximums
region_mins = regions(1,:);
region_maxs = regions(2,:);

shared_regions = regions; %We know shared regions size is limited by the regions size

%We want to loop through all points of interest both minimums and maximums,
%and label each by which pair they belong to. Once we find a point where
%all pairs have been "completed", that specifies a complete range

[sorted_points, inds] = sort([region_mins region_maxs]);
region_labels = ones(1,size(regions,2));
region_labels = [region_labels -1*region_labels]; %Positive for minimums and negative for maximums
region_labels = region_labels(inds);

M = 1; %M starts from 1
running_sum = 0;
for i = 1:numel(sorted_points)
    if running_sum == 0
        shared_regions(1,M) = sorted_points(i);
    end

    running_sum = running_sum + region_labels(i);

    if running_sum == 0
        shared_regions(2,M) = sorted_points(i);
        M = M + 1;
    end
end
%At the end all pairs should have been matched
assert(running_sum == 0);
M = M-1; %M was overincremented by 1

% M may have ended up smaller than N
shared_regions = shared_regions(:,1:M);

%There are edge cases where regions exactly touch, remove these:
shared_regions([false (shared_regions(2,1:(end-1)) == shared_regions(1,2:end))] ) = [];

end
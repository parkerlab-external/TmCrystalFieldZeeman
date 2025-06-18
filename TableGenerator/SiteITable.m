%Find the coarse lines given by energy levels

load(fullfile("..","Saves","coarse_levels.mat"));

gCFLabels = ["\pm 1/2" "\pm 3/2" "\pm 5/2" "\pm 7/2"];
eCFLabels = ["\pm 1/2" "\pm 3/2" "\pm 5/2"];

%Lines
levelLines = eE' - gE;
levelLineLabels = strcat("$", repmat(gCFLabels, [3 1]), "$", repmat(" to ", [3 4]), "$", repmat(eCFLabels', [1 4]), "$");

%Positions from amplitude fitting, and those used by simulation
load(fullfile("..","Saves","HFampfitting.mat"));
mFsource = round(mFsource, 2);
mFdest = round(mFdest, 2);

%Get centers
amp_fitted_centers = zeros(size(levelLines));
for i = 1:4:9
    for j = 1:4:13
        amp_fitted_centers((i-1)/4+1, (j-1)/4+1) = mean(mean(amp_fitted_positions(i:(i+3), j:(j+3))));
    end
end

%Hyperfine Lines
%Treat each 4x4 block separately
col1 = strings(0); %amp_fitted (observed)
col2 = strings(0); %level fitted
col3 = strings(0); %label (CF)
col4 = strings(0);
col5 = strings(0);
col6 = strings(0);
col7 = strings(0);

for i = 1:3
    for j = 1:4
        bigi = ((i-1)*4+1):((i-1)*4+4);
        bigj = ((j-1)*4+1):((j-1)*4+4);
        [hyp_positions, unique_inds] = uniquetol(amp_fitted_positions(bigi, bigj), 1e-12);

        header = true;
        

        for k = 1:numel(hyp_positions)
            hyp_pos = hyp_positions(k);
            degen_inds = abs(amp_fitted_positions - hyp_pos) < 1e-12;

            if k ==1
                col1(end+1) = sprintf("\\multirow{%d}{*}{%.4f}", numel(hyp_positions), amp_fitted_centers(i,j));
                col2(end+1) = sprintf("\\multirow{%d}{*}{%.4f}", numel(hyp_positions), levelLines(i,j));
                col3(end+1) = sprintf("\\multirow{%d}{*}{%s}", numel(hyp_positions), levelLineLabels(i,j));
            else
                col1(end+1) = "";
                col2(end+1) = "";
                col3(end+1) = "";
            end

            col4(end+1) = sprintf("%.4f", hyp_pos);
            col5(end+1) = sprintf("%.4f", sum(amp_fitted_amps(degen_inds), "all"));
            col6(end+1) = sprintf("%.4f", mean(amp_fitted_widths(degen_inds), "all"));

            %Labels
            hyp_source = unique(mFsource(degen_inds));
            if numel(hyp_source) == 2 && numel(unique(abs(hyp_source))) == 1
                hyp_source = sprintf("\\pm %d", unique(abs(hyp_source)));
            else
                hyp_source = join(string(hyp_source),",");
            end

            hyp_dest = unique(mFdest(degen_inds));
            if numel(hyp_dest) == 2 && numel(unique(abs(hyp_dest))) == 1
                hyp_dest = sprintf("\\pm %d", unique(abs(hyp_dest)));
            else
                hyp_dest = join(string(hyp_dest),",");
            end

            col7(end+1) = strcat("$", hyp_source, "$", " to ", "$",hyp_dest,"$");
        end
    end
end

fileID = fopen('SiteITable.txt','w');
BigTableString = [col1' col2' col3' col4' col5' col6' col7'];
for i = 1:size(BigTableString, 1)
    rowstring = join(BigTableString(i,:)," & ");
    fprintf(fileID, "%s", rowstring);
    fprintf(fileID, "\\\\\n");
end
fclose(fileID);

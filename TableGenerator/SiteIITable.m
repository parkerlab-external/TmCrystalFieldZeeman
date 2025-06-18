load(fullfile("..","Saves", "coarse_levels.mat"));
load(fullfile("..","Saves", "coarse_lines.mat"));

%For Site II manually make the labels
gCFLabels = flip(["\pm 1/2" "\pm 3/2" "\pm 5/2" "\pm 7/2"]);
eCFLabels = flip(["\pm 1/2" "\pm 3/2" "\pm 5/2"]);

%Lines
levelLines = eE' - gE;
levelLinesErr = sqrt(eE_err.^2' + gE_err.^2);
levelLineLabels = strcat("$", repmat(gCFLabels, [3 1]), "$", repmat(" to ", [3 4]), "$", repmat(eCFLabels', [1 4]), "$");
levelLines = levelLines';
levelLineLabels = levelLineLabels';
levelLinesErr = levelLinesErr';

fileID = fopen('SiteIITable.txt','w');
for i = 1:numel(levelLineLabels)
    my_line = levelLines(i);
    [~, line_ind] = min(abs(line_positions - my_line));

    level_pos = sprintf("%.4f",my_line);
    level_err = sprintf("%.4f",levelLinesErr(i));

    line_pos = sprintf("%.4f",line_positions(line_ind));
    line_width = sprintf("%.4f",line_widths(line_ind));
    line_err = sprintf("%.4f",line_positions_err(line_ind));

    if levelLineLabels(i) == "$\pm 7/2$ to $\pm 1/2$"
        line_pos = "";
        line_width = "";
        line_err = "";
    end
    
    fprintf(fileID, "%s", join([line_pos line_width line_err level_pos level_err levelLineLabels(i)], " & "));
    fprintf(fileID, "\\\\\n");
end
fclose(fileID)
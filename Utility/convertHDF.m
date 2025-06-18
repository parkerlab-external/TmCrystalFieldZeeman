FileList = dir(fullfile(cd, "**", "*.mat"));

for i = 1:numel(FileList)
    filepath = fullfile(FileList(i).folder, FileList(i).name);
    contents = load(filepath);
    save(filepath, "-struct", "contents", "-v7.3");
end
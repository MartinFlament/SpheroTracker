function list = listTiffs(folder)

if ~isdir(folder)
    error('Folder does not exist');
end
files=[];
files=dir(fullfile(folder,'*tif'));
if isempty(files)
files=dir(fullfile(folder,'*tiff'));
end
if isempty(files)
files=dir(fullfile(folder,'*TIF'));
end
if isempty(files)
files=dir(fullfile(folder,'*TIFF'));
end

list = {files.name}';

l = cellfun(@length, list);

if ~all(l(:) == l(1))
    warning('The filesnames have different lengths (trailing zeros?)');
end

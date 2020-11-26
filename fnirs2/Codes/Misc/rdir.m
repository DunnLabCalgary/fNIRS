function [files,folders,size] = rdir(absDir)
%Returns paths of all subfolders, files, subfolder files in a given
%directory, and size of the directory
%Input - (absDir) Absolute path to the directory of interest, ending in \
%Output - [files,folders,size] files: all files located in all absDir and
%subdirectories, folders: all folders located in all absDir and
%subdirectories, size: size in bytes of the absolute directory
fileList=dir(absDir);
folders={};
files={};
tempFiles={};
tempFolders={};
size=0;
if ~strcmp('\',absDir(end))
    absDir=[absDir '\'];
end
for i=1:length(fileList)
    if ~strcmp(fileList(i).name,'..') && ~strcmp(fileList(i).name,'.') && ~strncmp(fileList(i).name,'$',1)
        if isdir([absDir fileList(i).name])
            folders{end+1,1}=[absDir fileList(i).name '\'];
            if isempty(strfind([absDir fileList(i).name '\'],'.nff'))
                [tempFiles, tempFolders,tempSize]=rdir([absDir fileList(i).name '\']);
            end
            if ~isempty(tempFiles)
                files=[files; tempFiles];
                size=size+tempSize;
            end
            if ~isempty(tempFolders)
                folders=[folders; tempFolders];
            end
        else
            files{end+1,1}=[absDir fileList(i).name];
            size=size+fileList(i).bytes;
        end
    end
end
end
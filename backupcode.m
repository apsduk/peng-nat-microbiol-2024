function backupcode(fname, targetfolder)

% --- Determine file seperator --------------------------------------------
if ispc
    sep = '\';
elseif isunix
    sep = '/';
elseif ismac
    sep = '/';
else
    error('Unsupported OS');
end

% ----- Get list of dependencies ------------------------------------------
[fList, pList] = matlab.codetools.requiredFilesAndProducts([fname,'.m']);

% ----- Copy files to target folder and conver to txt ---------------------
for f = 1:length(fList)
    
    % ----- Split the string ----------------------------------------------
    if isempty(strfind(fList{f},'\'))
        txt = strsplit(fList{f}, '/');
    elseif isempty(strfind(fList{f},'/'))
        txt = strsplit(fList{f}, '\');
    end
    
    % ----- Select file name ----------------------------------------------
    filename = txt{end};
        
    % ----- Copy file -----------------------------------------------------    
    copyfile(fList{f}, [targetfolder,sep,filename(1:end-1),'txt']);
    
end

end
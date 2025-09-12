function copy_dependencies(src_dir,data_dir,funname)
% copies all .m file dependencies of the main.m file as .txt files as .txt
% files in the data_dir directory

% add .txt version of executed main.m file & dependent functions
fList = matlab.codetools.requiredFilesAndProducts(funname);

% exclude files from bin\external & files w/o .m extension
cd(src_dir); cd('..'); cd('bin'); cd('external'); binext_dir = pwd;
fList = fList(~startsWith(fList,binext_dir));
fList = fList(endsWith(fList,'.m'));

% make data\raw\subdir\mfiles directory into which to save used mfiles
cd(data_dir);
mkdir('mfiles');

% iterate over files
for kfL = 1:numel(fList)
    fn  = fList{kfL};
    [~,fn2,~] = fileparts(fn);
    fn2 = append(fn2,'.txt'); % .m -> .txt
    fn2 = fullfile(data_dir,'mfiles',fn2);
    copyfile(fn, fn2);
end

cd(src_dir);
end
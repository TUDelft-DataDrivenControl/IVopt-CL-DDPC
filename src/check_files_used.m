clc;
imp_files = {'plot_figs\plot_figs_paper.m','main.m','main_dRe.m','main_dp.m','main_dN.m','processing\main_processing.m'};
cd('src');
src_dir = pwd;
cd('..');
fList3 = {};
cd(src_dir);
for k = 1:length(imp_files)
fn = fullfile(pwd,imp_files{k});
fList = matlab.codetools.requiredFilesAndProducts(fn);
% exclude files from bin\external & files w/o .m extension
if k==1
    cd('..'); cd('bin'); bin_dir = pwd;
    cd(src_dir);
end
fList = fList(~startsWith(fList,bin_dir));
fList = fList(endsWith(fList,'.m'));
[~,fList2,~] = arrayfun(@(x) fileparts(x),fList);
fList3 = [fList3; fList2(:)];
clear fList2 fList
end
fList3 = unique(fList3);
%%
% for k = 1:length(fList3)
% fprintf('%s.m\n',fList3{k})
% end
for k = 1:length(fList3)
    fList3{k} = [fList3{k},'.m'];
end

%% all .m files
files = dir('**/*.m');
names = {files.name}';  % Column vector of filenames

%% unused .m files
unused = setdiff(names,fList3)
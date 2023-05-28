
clc
clear
close all

%%

directory  = [pwd, '\..\data\'];
cond       = {'tnt', 'sham'};

for i = 1:length(cond)
    
    dir1   = dir([directory, cond{i}]);
    n_sbj  = sum([dir1(~ismember({dir1.name},{'.','..'})).isdir]);

    for j = 1:n_sbj
        
        dir1   = [directory, cond{i}, '\Rat', num2str(j)];
        n_part = size(dir([dir1, '\*.txt']), 1);
        
        data = [];
        for k = 1:n_part
            tmp = load([dir1, '\', num2str(k), '.txt']);
            data = [data; tmp];
        end
        
        save([dir1, '\Rat', num2str(j)], 'data');
    end
    
end


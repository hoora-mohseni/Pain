clc
clear
close 

%%

n = 30;
t = 2*n

%%

load('ttl');
load('ssl');


idx = randperm(size(ttl, 2));
ttl = ttl(:, idx);
ttl = ttl(:, 1:size(ssl, 2));

idx = randperm(size(ttl, 2));
idx = idx(1:n/2);

ttlu = ttl(:, idx);
sslu = ssl(:, idx);

data  = [ttlu, sslu];
label = [ones(1, n/2), zeros(1, n/2)]';

idx   = randperm(n);

label = label(idx);
data  = data(:, idx);
data  = data(:);

fileID = fopen(['data_',num2str(t),'s.txt'],'w');
fileID2 = fopen(['label_',num2str(t),'s.txt'],'w');

for i = 1:length(data)
    fprintf(fileID, [num2str(abs(data(i))), '\n']);
end

for i = 1:length(label)
    fprintf(fileID2, [num2str(label(i)), '\n']);
end





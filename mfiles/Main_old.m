
clc
clear
close all

%% loading dataset

directory  = [pwd, '\..\data\'];
cond       = {'tnt', 'sham'};
n_sbj      = 5;

tnt  = cell(1, n_sbj);
sham = cell(1, n_sbj);
for i = 1:2 % number of conditions (tnt, sham)
    dir1   = dir([directory, cond{i}]);
    for j = 1:n_sbj % number of subjects in each condition
        name        = [directory, cond{i}, '\Rat', num2str(j), '\Rat', num2str(j)];
        if i == 1
           tnt{j}  = load(name).data; 
        else
           sham{j} = load(name).data; 
        end
    end % j
end %i

tntf   = tnt;
shamf  = sham;

%% filtering

fs           = 512;
cutoff_low   = 0.5;
cutoff_high  = 60;
filter_order = 6;

%-------- low pass filter ----------%
Wn           = cutoff_high/(fs/2);             % Normalized cutoff frequency
[Bh, Ah]     = butter(filter_order, Wn, 'low');           % Butterworth

%-------- high pass filter ----------%
Wn           = cutoff_low/(fs/2);              % Normalized cutoff frequency
[Bl, Al]     = butter(filter_order, Wn, 'high');          % Butterworth

%-------- filtering ----------%
for i = 1:n_sbj
    
    % lowpass
    tntf{i}(:, 2:end)  = filtfilt(Bh, Ah, tnt{i}(:, 2:end));
    shamf{i}(:, 2:end) = filtfilt(Bh, Ah, sham{i}(:, 2:end));
    
    % highpass
    tntf{i}(:, 2:end)  = filtfilt(Bl, Al, tntf{i}(:, 2:end));
    shamf{i}(:, 2:end) = filtfilt(Bl, Al, shamf{i}(:, 2:end));
    
end


%% Validation

% a = fftshift(abs(fft(tnt{1}(:, 2))));
% b = fftshift(abs(fft(tntf{1}(:, 2))));
% f = -fs/2:fs/length(tnt{1}(:, 2)):fs/2-fs/length(tnt{1}(:, 2));
% 
% plot(f, a); hold on; plot(f, b)



clc
clear
close all

%% Loading data

dirs    = dir([pwd, '/../data/*/*/*.txt']);
n_files = length(dirs);


for i = 1:n_files
    if i < 10
        fprintf([num2str(i), '    | ', dirs(i).name, '\n']);
    else
        fprintf([num2str(i), '   | ', dirs(i).name, '\n']);
    end
end

inp = input('\n\nPlease enter the number associated with your file: ');

while(inp > n_files || inp < 1 || inp ~= round(inp))
    inp = input('Please TRY AGAIN!!! Enter the number associated with your file: ');
end

data = load([dirs(inp).folder, '/', dirs(inp).name]);
name = dirs(inp).name;

ts_in = input('start time (s) = ');
t_in  = input('length of signal (s) = ');

fs           = 512;
data = data(round(ts_in*fs) + (1:round(fs*t_in)), :);

%% filtering


cutoff_low   = 0.5;
cutoff_high  = 60;
filter_order = 6;

%-------- low pass filter ----------%
Wn           = cutoff_high/(fs/2);             % Normalized cutoff frequency
[Bh, Ah]     = butter(filter_order, Wn, 'low');           % Butterworth

%-------- high pass filter ----------%
Wn           = cutoff_low/(fs/2);              % Normalized cutoff frequency
[Bl, Al]     = butter(filter_order, Wn, 'high');          % Butterworth

dataf = filtfilt(Bl, Al, data);
dataf = filtfilt(Bh, Ah, dataf);

%% seperating band 


cutoff3      = 3.5; 
cutoff7      = 7;
cutoff14     = 14;
cutoff30     = 30;

%-------- delta band ----------%

Wn           = cutoff3/(fs/2);            
[B3, A3]     = butter(filter_order, Wn, 'low');  

dataf_delta = filtfilt(B3, A3, dataf);

%-------- theta band ----------%


Wn           = cutoff3/(fs/2);            
[B3, A3]     = butter(filter_order, Wn, 'high');  

Wn           = cutoff7/(fs/2); 
[B7, A7]     = butter(filter_order, Wn, 'low'); 

dataf_theta  = filtfilt(B3, A3, dataf);
dataf_theta  = filtfilt(B7, A7, dataf_theta);

%-------- alfa band ----------%

Wn           = cutoff7/(fs/2);            
[B7, A7]     = butter(filter_order, Wn, 'high');  

Wn           = cutoff14/(fs/2); 
[B14, A14]   = butter(filter_order, Wn, 'low'); 


dataf_alpha  = filtfilt(B7, A7, dataf);
dataf_alpha  = filtfilt(B14, A14, dataf_alpha);


%-------- beta band ----------%


Wn           = cutoff14/(fs/2);            
[B14, A14]   = butter(filter_order, Wn, 'high');  

Wn           = cutoff30/(fs/2); 
[B30, A30]   = butter(filter_order, Wn, 'low'); 

dataf_betha  = filtfilt(B14, A14, dataf);
dataf_betha  = filtfilt(B30, A30, dataf_betha);


%-------- gamma band ----------%

Wn           = cutoff30/(fs/2);            
[B30, A30]   = butter(filter_order, Wn, 'high');  

dataf_gamma  = filtfilt(B30, A30, dataf);

%%

fs = 512;
time = data(:, 1);

%%

for b = 1:7
    
    aaa = 'signal-raw';
    
    switch b

        case 2
            data = dataf;
            aaa = 'signal-filtered';
        case 3
            data = dataf_delta;
            aaa = 'delta-band';
        case 4
            data = dataf_theta;
            aaa = 'theta band';
        case 5
            data = dataf_alpha;
            aaa = 'alpha band';
        case 6
            data = dataf_betha;
            aaa = 'beta-band';
        case 7
            data = dataf_gamma;
            aaa = 'gamma-band';
    end
    
    data(:, 1) = time;
    
    data_fft = abs(fft(data));
    fff      = -fs/2:fs/size(data, 1):fs/2-fs/size(data, 1);

    figure('color', 'w');
    for i = 2:size(data, 2)
       subplot(2, 2, i-1)
       periodogram(data(:, i),rectwin(size(data, 1)), size(data, 1), fs);
       title(['Channel ', num2str(i-1), ' - ', aaa]);
       box on, ylabel('Power/Freq (db/Hz)'); xlabel('Frequency (Hz)');
    end


    figure('color', 'w');
    for i = 2:size(data, 2)
       subplot(2, 2, i-1)
       plot(data(:, 1)/60, data(:, i));
       title(['Channel ', num2str(i-1), ' - ', aaa]);
       box on, ylabel('Amplitude'); xlabel('Time (Min)'); grid on
    end

end

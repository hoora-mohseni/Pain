
clc
clear
close all

%% loading dataset

directory      = [pwd, '\..\data\'];
cond           = {'tnt', 'sham'};
n_sbj          = 5;
n_channel      = 4;
lowerTime      = 0;     % seconds
upperTime      = 30;    % seconds

tnt  = cell(1, n_sbj);
sham = cell(1, n_sbj);
for i = 1:2 % number of conditions (tnt, sham)
    dir1   = dir([directory, cond{i}]);
    for j = 1:n_sbj % number of subjects in each condition
        name        = [directory, cond{i}, '\Rat', num2str(j), '\Rat', num2str(j)];
        if i == 1
            t       = load(name);
            tnt{j}  = t.data    ; 
            
        else
           sh       = load(name); 
           sham{j}  = sh.data   ;
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


%% seperating band 

% cut off freq.
cutoff3      = 3.5; 
cutoff7      = 7;
cutoff14     = 14;
cutoff30     = 30;

%-------- delta band ----------%

Wn           = cutoff3/(fs/2);            
[B3, A3]     = butter(filter_order, Wn, 'low');           

tntdelta     = tntf;
shamdelta    = shamf;
for i = 1:n_sbj
    tntdelta{i}(:, 2:end)  = filtfilt(B3, A3, tntf{i}(:, 2:end));
    shamdelta{i}(:, 2:end) = filtfilt(B3, A3, shamf{i}(:, 2:end));
end 

%-------- theta band ----------%

tnttheta  = tntf ;
shamtheta = shamf;

Wn           = cutoff3/(fs/2);            
[B3, A3]     = butter(filter_order, Wn, 'high');  

Wn           = cutoff7/(fs/2); 
[B7, A7]     = butter(filter_order, Wn, 'low'); 

for i = 1:n_sbj
    
    tnttheta{i}(:, 2:end)  = filtfilt(B3, A3, tntf{i}(:, 2:end));
    shamtheta{i}(:, 2:end) = filtfilt(B3, A3, shamf{i}(:, 2:end));
    
    
    tnttheta{i}(:, 2:end)  = filtfilt(B7, A7, tnttheta{i}(:, 2:end));
    shamtheta{i}(:, 2:end) = filtfilt(B7, A7, shamtheta{i}(:, 2:end));
    
    
end 

%-------- alpha band ----------%

tntalpha  = tntf ;
shamalpha = shamf;

Wn           = cutoff7/(fs/2);            
[B7, A7]     = butter(filter_order, Wn, 'high');  

Wn           = cutoff14/(fs/2); 
[B14, A14]   = butter(filter_order, Wn, 'low'); 

for i = 1:n_sbj
    
    tntalpha{i}(:, 2:end)  = filtfilt(B7, A7, tntf{i}(:, 2:end));
    shamalpha{i}(:, 2:end) = filtfilt(B7, A7, shamf{i}(:, 2:end));
    
    
    tntalpha{i}(:, 2:end)  = filtfilt(B14, A14, tntalpha{i}(:, 2:end));
    shamalpha{i}(:, 2:end) = filtfilt(B14, A14, shamalpha{i}(:, 2:end));
    
    
end 

%-------- beta band ----------%

tntbeta  = tntf ;
shambeta = shamf;

Wn           = cutoff14/(fs/2);            
[B14, A14]     = butter(filter_order, Wn, 'high');  

Wn           = cutoff30/(fs/2); 
[B30, A30]     = butter(filter_order, Wn, 'low'); 

for i = 1:n_sbj
    
    tntbeta{i}(:, 2:end)  = filtfilt(B14, A14, tntf{i}(:, 2:end));
    shambeta{i}(:, 2:end) = filtfilt(B14, A14, shamf{i}(:, 2:end));
    
    
    tntbeta{i}(:, 2:end)  = filtfilt(B30, A30, tntbeta{i}(:, 2:end));
    shambeta{i}(:, 2:end) = filtfilt(B30, A30, shambeta{i}(:, 2:end));
    
    
end 

%-------- gamma band ----------%

tntgamma  = tntf ;
shamgamma = shamf;

Wn           = cutoff30/(fs/2);            
[B30, A30]     = butter(filter_order, Wn, 'high');  


for i = 1:n_sbj
    
    tntgamma{i}(:, 2:end)  = filtfilt(B30, A30, tntf{i}(:, 2:end));
    shamgamma{i}(:, 2:end) = filtfilt(B30, A30, shamf{i}(:, 2:end));
    
   
end 

%% ----figure----- 

for rep = 1:2
    
    if rep == 1
        show      = shamf;
        showalpha = shamalpha;
        showbeta  = shambeta;
        showdelta = shamdelta;
        showgamma = shamgamma;
        showtheta = shamtheta;
    else
        show      = tntf;
        showalpha = tntalpha;
        showbeta  = tntbeta;
        showdelta = tntdelta;
        showgamma = tntgamma;
        showtheta = tnttheta;
    end
    
    for i = 1:n_sbj
        for j = 2:n_channel+1

            f1 = figure('WindowState','maximized'); 
            sgtitle([upper(cond{rep}), ' - Rat ', num2str(i), ' - Channel ', num2str(j-1)]);
            
            subplot(3,2,1);
            plot(show{i}(:,1), show{i}(:,j));
%             xlabel('time (seconds)');
            ylabel('Amplitude (mV)');
            xlim([lowerTime, upperTime]);
            title('Channel');
            box off
            
            
            subplot(3,2,4);
            plot(show{i}(:,1), showalpha{i}(:,j));
%             xlabel('time (seconds)');
            ylabel('Amplitude (mV)');
            xlim([lowerTime, upperTime]);
            box off
            title('Alpha');
            
            subplot(3,2,5);
            plot(show{i}(:,1), showbeta{i}(:,j));
            xlabel('time (seconds)');
            ylabel('Amplitude (mV)');
            xlim([lowerTime, upperTime]);
            box off
            title('Beta');
            
            subplot(3,2,2);
            plot(show{i}(:,1), showdelta{i}(:,j));
%             xlabel('time (seconds)');
            ylabel('Amplitude (mV)');
            xlim([lowerTime, upperTime]);
            box off
            title('Delta');
            
            subplot(3,2,3);
            plot(show{i}(:,1), showtheta{i}(:,j));
%             xlabel('time (seconds)');
            ylabel('Amplitude (mV)');
            box off
            xlim([lowerTime, upperTime]);
            title('Theta');
            
            subplot(3,2,6);
            plot(show{i}(:,1), showgamma{i}(:,j));
            xlabel('time (seconds)');
            ylabel('Amplitude (mV)');
            box off
            xlim([lowerTime, upperTime]);
            title('Gamma');
            
            set(gcf, 'color', 'w');
            
            pause(3);
            f1.Position;
            
            saveas(gcf, [directory, cond{rep}, '\Rat', num2str(i), '\Rat', num2str(i), '_Ch', num2str(j-1), '_', num2str(lowerTime), 'to', num2str(upperTime),'.jpeg']);
            
            close all
        end
    end
end


% a = fftshift(abs(fft(tntf{1}(:, 2))));
% b = fftshift(abs(fft(tntgamma{1}(:, 2))));
% f = -fs/2:fs/length(tntf{1}(:, 2)):fs/2-fs/length(tntf{1}(:, 2));
% 
% plot(f, a); hold on; plot(f, b)

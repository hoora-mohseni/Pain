
clc
clear
close all


%% Constants

wl        = 3;
tint      = 30;

IMF_n     = 2;

%% loading dataset

directory  = [pwd, '/../data/'];
cond       = {'tnt', 'sham'};
n_sbj      = 5;

tnt  = cell(1, n_sbj);
sham = cell(1, n_sbj);
for i = 1:2 % number of conditions (tnt, sham)
    dir1   = dir([directory, cond{i}]);
    for j = 1:n_sbj % number of subjects in each condition
        name        = [directory, cond{i}, '/Rat', num2str(j), '/Rat', num2str(j)];
        if i == 1
            t       = load(name);
            tnt{j}  = t.data(:, 2:end); 
            
        else
           sh       = load(name); 
           sham{j}  = sh.data(:, 2:end);
        end
    end % j
end %i

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
% % % for i = 1:n_sbj
% % %     
% % %     % lowpass
% % %     tnt{i}  = filtfilt(Bh, Ah, tnt{i});
% % %     sham{i} = filtfilt(Bh, Ah, sham{i});
% % %     
% % %     % highpass
% % %     tnt{i}  = filtfilt(Bl, Al, tnt{i});
% % %     sham{i} = filtfilt(Bl, Al, sham{i});
% % %     
% % % end


%% Validation

% a = fftshift(abs(fft(tnt{1}(:, 2))));
% b = fftshift(abs(fft(tntf{1}(:, 2))));
% f = -fs/2:fs/length(tnt{1}(:, 2)):fs/2-fs/length(tnt{1}(:, 2));
% 
% plot(f, a); hold on; plot(f, b)


%% seperating band 


cutoff3      = 3.5; 
cutoff7      = 7;
cutoff14     = 14;
cutoff30     = 30;

%-------- delta band ----------%

tntdelta     = tnt;
shamdelta    = sham;

Wn           = cutoff3/(fs/2);            
[B3, A3]     = butter(filter_order, Wn, 'low');           

for i = 1:n_sbj  
    tntdelta{i}  = filtfilt(B3, A3, tnt{i});
    shamdelta{i} = filtfilt(B3, A3, sham{i});   
end 

%-------- theta band ----------%

tnttheta  = tnt ;
shamtheta = sham;

Wn           = cutoff3/(fs/2);            
[B3, A3]     = butter(filter_order, Wn, 'high');  

Wn           = cutoff7/(fs/2); 
[B7, A7]     = butter(filter_order, Wn, 'low'); 

for i = 1:n_sbj 
    tnttheta{i}  = filtfilt(B3, A3, tnt{i});
    shamtheta{i} = filtfilt(B3, A3, sham{i});
        
    tnttheta{i}  = filtfilt(B7, A7, tnttheta{i});
    shamtheta{i} = filtfilt(B7, A7, shamtheta{i});   
end 

%-------- alfa band ----------%

tntalfa  = tnt ;
shamalfa = sham;

Wn           = cutoff7/(fs/2);            
[B7, A7]     = butter(filter_order, Wn, 'high');  

Wn           = cutoff14/(fs/2); 
[B14, A14]   = butter(filter_order, Wn, 'low'); 

for i = 1:n_sbj
    tntalfa{i}  = filtfilt(B7, A7, tnt{i});
    shamalfa{i} = filtfilt(B7, A7, sham{i});
    
    tntalfa{i}  = filtfilt(B14, A14, tntalfa{i});
    shamalfa{i} = filtfilt(B14, A14, shamalfa{i});
end 

%-------- beta band ----------%

tntbeta  = tnt ;
shambeta = sham;

Wn           = cutoff14/(fs/2);            
[B14, A14]   = butter(filter_order, Wn, 'high');  

Wn           = cutoff30/(fs/2); 
[B30, A30]   = butter(filter_order, Wn, 'low'); 

for i = 1:n_sbj
    tntbeta{i}  = filtfilt(B14, A14, tnt{i});
    shambeta{i} = filtfilt(B14, A14, sham{i});
    
    tntbeta{i}  = filtfilt(B30, A30, tntbeta{i});
    shambeta{i} = filtfilt(B30, A30, shambeta{i}); 
end 

%-------- gamma band ----------%

tntgamma  = tnt ;
shamgamma = sham;

Wn           = cutoff30/(fs/2);            
[B30, A30]   = butter(filter_order, Wn, 'high');  


for i = 1:n_sbj
    tntgamma{i}  = filtfilt(B30, A30, tnt{i});
    shamgamma{i} = filtfilt(B30, A30, sham{i});
end 


TNT  = cell(size(tnt));
SHAM = cell(size(sham));
for i = 1:numel(tnt)
    TNT{i}(:, :, 1)  = tnt{i};
    TNT{i}(:, :, 6)  = tntdelta{i};
    TNT{i}(:, :, 5)  = tnttheta{i};
    TNT{i}(:, :, 4)  = tntalfa{i};
    TNT{i}(:, :, 3)  = tntbeta{i};
    TNT{i}(:, :, 2)  = tntgamma{i};
    
    SHAM{i}(:, :, 1) = sham{i};
    SHAM{i}(:, :, 6) = shamdelta{i};
    SHAM{i}(:, :, 5) = shamtheta{i};
    SHAM{i}(:, :, 4) = shamalfa{i};
    SHAM{i}(:, :, 3) = shambeta{i};
    SHAM{i}(:, :, 2) = shamgamma{i};
end

%% -------seperation data-------

tntltime  = [106, 65 , 73 , 290, 71  ; 
             172, 187, 219, 356, 140 ;
             245, 260, 291, 428, 211];
         
tntrtime  = [440, 379, 365, 522, 299 ; 
             511, 455, 427, 597, 370 ;
             578, 524, 510, 668, 463];
         
shamltime = [68 , 87 , 161, 65 , 64  ;
             145, 152, 214, 126, 128 ; 
             210, 215, nan, 195, 194];
         
shamrtime = [277, 290, nan, 265, 272 ;
             460, 362, nan, 324, 335 ;
             525, 447, nan, 385, 401];


for i = 1:n_sbj
    cnt  = 0;
    cnt1 = 0;
    cnt2 = 0;
    for dose = 1:3
        for tim = (0:wl:tint-wl)*fs
             
             cnt = cnt + 1;
             TNTleft{i}(:, cnt, :, :)   = TNT{i}(round(tntltime(dose, i)*fs+tim+(0:wl*fs-1)), :, :);
             TNTright{i}(:, cnt, :, :)  = TNT{i}(round(tntrtime(dose, i)*fs+tim+(0:wl*fs-1)), :, :);
             
             if isnan(shamltime(dose, i))
                 a = 0;
             end
             
             if ~isnan(shamltime(dose, i))
                 cnt1 = cnt1 + 1;
                 SHAMleft{i}(:, cnt1, :, :)  = SHAM{i}(round(shamltime(dose, i)*fs+tim+(0:wl*fs-1)), :, :);
             end
             
             if ~isnan(shamrtime(dose, i))
                 cnt2 = cnt2 + 1;
                 SHAMright{i}(:, cnt2, :, :) = SHAM{i}(round(shamrtime(dose, i)*fs+tim+(0:wl*fs-1)), :, :);
             end
             
        end
    end
end

%-----------


TL = TNTleft{1};
for i = 2:numel(TNTleft)
    [~, N, ~, ~] = size(TNTleft{i});
    TL(:, end+1:end+N, :, :)  = TNTleft{i};
end
[M, N, O, P] = size(TL);
TL    = reshape(TL, M, N, O*P);
TLimf = TL;

for i = 1:N
    for j = 1:O*P-8
        IMF = emd(TLimf(:, i, j));
        TLimf(:, i, j, 1+(1:IMF_n)) = IMF(:, 1:IMF_n);
    end
end
[M, N, O, P] = size(TLimf);
TLimf       = reshape(TLimf, M, N, O*P);



%-----------


SL = SHAMleft{1};
for i = 2:numel(SHAMleft)
    [~, N, ~, ~] = size(SHAMleft{i});
    SL(:, end+1:end+N, :, :)  = SHAMleft{i};
end
[M, N, O, P] = size(SL);
SL  = reshape(SL, M, N, O*P);
SLimf = SL;

for i = 1:N
    for j = 1:O*P-8
        IMF = emd(SLimf(:, i, j));
        SLimf(:, i, j, 1+(1:IMF_n)) = IMF(:, 1:IMF_n);
    end
end
[M, N, O, P] = size(SLimf);
SLimf       = reshape(SLimf, M, N, O*P);

%%

shamm = cell(1, 5);
for i = 1:5
    for j = 1:size(SHAMleft{i}, 2)
        kk = 0;
        for k = [1, 3]
            kk = kk+1;
            tmp = squeeze(SHAMleft{i}(:, j, kk, 1));
            shamm{i}((j-1)*numel(tmp)+1:j*numel(tmp), kk) = tmp;
        end
    end
end


tntt = cell(1, 5);
for i = 1:5
    for j = 1:size(TNTleft{i}, 2)
        kk = 0;
        for k = [1, 3]
            kk  = kk+1;
            tmp = squeeze(TNTleft{i}(:, j, kk, 1));
            tntt{i}((j-1)*numel(tmp)+1:j*numel(tmp), kk) = tmp;
        end
    end
end

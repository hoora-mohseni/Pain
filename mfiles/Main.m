
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
for i = 1:n_sbj
    
    % lowpass
    tnt{i}  = filtfilt(Bh, Ah, tnt{i});
    sham{i} = filtfilt(Bh, Ah, sham{i});
    
    % highpass
    tnt{i}  = filtfilt(Bl, Al, tnt{i});
    sham{i} = filtfilt(Bl, Al, sham{i});
    
end



%% Validation

% a = fftshift(abs(fft(tnt{1}(:, 2))));
% b = fftshift(abs(fft(tntf{1}(:, 2))));
% f = -fs/2:fs/length(tnt{1}(:, 2)):fs/2-fs/length(tnt{1}(:, 2));
% 
% plot(f, a); hold on; plot(f, b)


%% seperating band 

% for i = 1:5
%     [xxx, yyy] = find(tnt{i} >=1);
%     tnt{i}(xxx, yyy) = 0.99;
% 
%     [xxx, yyy] = find(sham{i} >=1);
%     sham{i}(xxx, yyy) = 0.99;
%     
%     [xxx, yyy] = find(tnt{i} <=-1);
%     tnt{i}(xxx, yyy) = -0.99;
% 
%     [xxx, yyy] = find(sham{i} <=-1);
%     sham{i}(xxx, yyy) = -0.99;
%     
% end


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


% 
% for i = 1:5
%     
%     for j = 1:5
%         
%         switch j 
%             
%             case 1
%                 ttttt = tntdelta;
%                 shsh  = shamdelta;
%             case 2
%                 ttttt = tnttheta;
%                 shsh  = shamtheta;
%             case 3
%                 ttttt = tntalfa;
%                 shsh  = shamalfa;
%             case 4
%                 ttttt = tntbeta;
%                 shsh  = shambeta;
%             case 5
%                 ttttt = tntgamma;
%                 shsh  = shamgamma;
%         end
%         
%     
%         [xxx, yyy] = find(ttttt{i} >=1);
%         ttttt{i}(xxx, yyy) = 0.9999999;
%        
%         [xxx, yyy] = find(ttttt{i} <=-1);
%         ttttt{i}(xxx, yyy) = -0.9999999;
%         
%         [xxx, yyy] = find(shsh{i} >=1);
%         shsh{i}(xxx, yyy) = 0.9999999;
% 
%         [xxx, yyy] = find(shsh{i} <=-1);
%         shsh{i}(xxx, yyy) = -0.9999999;
%         
%         
%         switch j 
%             
%             case 1
%                 tntdelta = ttttt;
%                 shamdelta = shsh;
%             case 2
%                 tnttheta = ttttt;
%                 shamtheta = shsh;
%             case 3
%                 tntalfa = ttttt;
%                 shamalfa = shsh;
%             case 4
%                 tntbeta = ttttt;
%                 shambeta = shsh;
%             case 5
%                 tntgamma = ttttt;
%                 shamgamma = shsh;
%         end
%         
%         
% 
%     end
% end


% for i = 1:5
%     tnt{i}      = round(10e4*tnt{i});
%     tntdelta{i} = round(10e4*tntdelta{i});
%     tnttheta{i} = round(10e4*tnttheta{i});
%     tntalfa{i}  = round(10e4*tntalfa{i});
%     tntbeta{i}  = round(10e4*tntbeta{i});
%     tntgamma{i} = round(10e4*tntgamma{i});
%     
%     sham{i}      = round(10e4*sham{i});
%     shamdelta{i} = round(10e4*shamdelta{i});
%     shamtheta{i} = round(10e4*shamtheta{i});
%     shamalfa{i}  = round(10e4*shamalfa{i});
%     shambeta{i}  = round(10e4*shambeta{i});
%     shamgamma{i} = round(10e4*shamgamma{i});
%     
% end


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


TR = TNTright{1};
for i = 2:numel(TNTright)
    [~, N, ~, ~] = size(TNTright{i});
    TR(:, end+1:end+N, :, :)  = TNTright{i};
end

[M, N, O, P] = size(TR);
TR    = reshape(TR, M, N, O*P);
TRimf = TR;

for i = 1:N
    for j = 1:O*P-8
        IMF = emd(TRimf(:, i, j));
        TRimf(:, i, j, 1+(1:IMF_n)) = IMF(:, 1:IMF_n);
    end
end
[M, N, O, P] = size(TRimf);
TRimf    = reshape(TRimf, M, N, O*P);


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

%-----------


SR = SHAMright{1};
for i = 2:numel(SHAMright)
    [~, N, ~, ~] = size(SHAMright{i});
    SR(:, end+1:end+N, :, :)  = SHAMright{i};
end
[M, N, O, P] = size(SR);
SR  = reshape(SR, M, N, O*P);
SRimf = SR;

for i = 1:N
    for j = 1:O*P-8
        IMF = emd(SRimf(:, i, j));
        SRimf(:, i, j, 1+(1:IMF_n)) = IMF(:, 1:IMF_n);
    end
end
[M, N, O, P] = size(SRimf);
SRimf    = reshape(SRimf, M, N, O*P);


%% Features F

for i = 1:size(SL, 2)
    for j = 1:size(SL, 3)
        shannon_entropy(i, j)   = wentropy(SL(:, i, j), 'shannon');
        mean_frequency(i, j)    = meanfreq(SL(:, i, j), fs);
        spectral_entropy(i, j)  = sum(pentropy(SL(:, i, j), fs));
        tmp                     = aryule(SL(:, i, j), 3);
        AR_coeff(i, j, 1:3)     = tmp(2:4);
    end
end
[M, N, O] = size(AR_coeff);
AR_coeff  = reshape(AR_coeff, M, O*N);

F_SL = [squeeze(mean(SL)), ... % mean              
        squeeze(var(SL)), ...  % variance
        squeeze(mode(SL)), ... % mode
        squeeze(median(SL)),...   % median
        squeeze(kurtosis(SL)),... % kurotsis
        squeeze(skewness(SL)),... % skewness 
        squeeze(std(SL)),...      %standard deviation
        squeeze(std(SL))./squeeze(mean(SL)), ... % coefficient of variation
        squeeze(quantile(SL, .25)),... % first quartile
        squeeze(quantile(SL, .75)),... % third quartile
        squeeze(iqr(SL)),...    % interquartile range
        squeeze(rms(SL)),...    % root mean square
        squeeze(rms(SL).^2),... % average power
        squeeze(sum(SL.^2)),... % energy
        shannon_entropy, ...    % shannon entropy
        mean_frequency,...      % mean frequency
        spectral_entropy,...    % spectral entropy
        AR_coeff];              % AR coefficients           


clear shannon_entropy mean_frequency spectral_entropy tmp AR_coeff



for i = 1:size(SR, 2)
    for j = 1:size(SR, 3)
        shannon_entropy(i, j)   = wentropy(SR(:, i, j), 'shannon');
        mean_frequency(i, j)    = meanfreq(SR(:, i, j), fs);
        spectral_entropy(i, j)  = sum(pentropy(SR(:, i, j), fs));
        tmp                     = aryule(SR(:, i, j), 3);
        AR_coeff(i, j, 1:3)     = tmp(2:4);
    end
end
[M, N, O] = size(AR_coeff);
AR_coeff    = reshape(AR_coeff, M, O*N);

F_SR = [squeeze(mean(SR)), ... % mean              
        squeeze(var(SR)), ...  % variance
        squeeze(mode(SR)), ... % mode
        squeeze(median(SR)),...   % median
        squeeze(kurtosis(SR)),... % kurotsis
        squeeze(skewness(SR)),... % skewness 
        squeeze(std(SR)),...    % standard deviation
        squeeze(std(SR))./squeeze(mean(SR)), ... % coefficient of variation
        squeeze(quantile(SR, .25)),... % first quartile
        squeeze(quantile(SR, .75)),... % third quartile
        squeeze(iqr(SR)),...    % interquartile range
        squeeze(rms(SR)),...    % root mean square
        squeeze(rms(SR).^2),... % average power
        squeeze(sum(SR.^2)),... % energy
        shannon_entropy, ...    % shannon entropy
        mean_frequency,...      % mean frequency
        spectral_entropy,...    % spectral entropy
        AR_coeff];              % AR coefficients   

clear shannon_entropy mean_frequency spectral_entropy tmp AR_coeff

for i = 1:size(TR, 2)
    for j = 1:size(TR, 3)
        shannon_entropy(i, j)   = wentropy(TR(:, i, j), 'shannon');
        mean_frequency(i, j)    = meanfreq(TR(:, i, j), fs);
        spectral_entropy(i, j)  = sum(pentropy(TR(:, i, j), fs));
        tmp                     = aryule(TR(:, i, j), 3);
        AR_coeff(i, j, 1:3)     = tmp(2:4);
    end
end
[M, N, O] = size(AR_coeff);
AR_coeff    = reshape(AR_coeff, M, O*N);

F_TR = [squeeze(mean(TR)), ... % mean              
        squeeze(var(TR)), ...  % variance
        squeeze(mode(TR)), ... % mode
        squeeze(median(TR)),...   % median
        squeeze(kurtosis(TR)),... % kurotsis
        squeeze(skewness(TR)),... % skewness 
        squeeze(std(TR)),...    % standard deviation
        squeeze(std(TR))./squeeze(mean(TR)), ... % coefficient of variation
        squeeze(quantile(TR, .25)),... % first quartile
        squeeze(quantile(TR, .75)),... % third quartile
        squeeze(iqr(TR)),...    % interquartile range
        squeeze(rms(TR)),...    % root mean square
        squeeze(rms(TR).^2),... % average power
        squeeze(sum(TR.^2)),... % energy
        shannon_entropy, ...    % shannon entropy
        mean_frequency,...      % mean frequency
        spectral_entropy,...    % spectral entropy
        AR_coeff];              % AR coefficients   


clear shannon_entropy mean_frequency spectral_entropy tmp AR_coeff

for i = 1:size(TL, 2)
    for j = 1:size(TL, 3)
        shannon_entropy(i, j)   = wentropy(TL(:, i, j), 'shannon');
        mean_frequency(i, j)    = meanfreq(TL(:, i, j), fs);
        spectral_entropy(i, j)  = sum(pentropy(TL(:, i, j), fs));
        tmp                     = aryule(TL(:, i, j), 3);
        AR_coeff(i, j, 1:3)     = tmp(2:4);
    end
end
[M, N, O] = size(AR_coeff);
AR_coeff  = reshape(AR_coeff, M, O*N);

F_TL = [squeeze(mean(TL)), ... % mean              
        squeeze(var(TL)), ...  % variance
        squeeze(mode(TL)), ... % mode
        squeeze(median(TL)),...   % median
        squeeze(kurtosis(TL)),... % kurotsis
        squeeze(skewness(TL)),... % skewness 
        squeeze(std(TL)),...    % standard deviation
        squeeze(std(TL))./squeeze(mean(TL)), ... % coefficient of variation
        squeeze(quantile(TL, .25)),... % first quartile
        squeeze(quantile(TL, .75)),... % third quartile
        squeeze(iqr(TL)),...    % interquartile range
        squeeze(rms(TL)),...    % root mean square
        squeeze(rms(TL).^2),... % average power
        squeeze(sum(TL.^2)),... % energy
        shannon_entropy, ...    % shannon entropy
        mean_frequency,...      % mean frequency
        spectral_entropy,...    % spectral entropy
        AR_coeff];              % AR coefficients   


clear shannon_entropy mean_frequency spectral_entropy tmp AR_coeff


F_SL(:, end+1) =   zeros(size(F_SL, 1), 1);
F_SR(:, end+1) =   ones(size(F_SR, 1), 1);
F_TL(:, end+1) = 2*ones(size(F_TL, 1), 1);
F_TR(:, end+1) = 3*ones(size(F_TR, 1), 1);

F    = [F_SL; F_SR; F_TL; F_TR];


band = {'total', 'gamma', 'betha', 'alpha', 'theta', 'delta'};
chan = {'1', '2', '3', '4'};
feat = {'mean', 'variance', 'mode', 'median', 'kurtosis', ...
        'skewness', 'SD', 'COV', 'fisrtQ', 'thirdQ', 'IQR', 'RMS', ...
        'POWER', 'ENERGY', 'ShanonEntropy', 'meanFreq', 'SpectralEntropy', 'AR1', 'AR2', 'AR3'};

kk = 0;
for i = 1:length(feat)
    for j = 1:length(band)
        for k = 1:length(chan)
            kk = kk +1;
            LABELL{kk} = [feat{i}, '-ch', chan{k}, '-', band{j}];
        end
    end
end

save([pwd, '/../F_', num2str(tint), '_', num2str(wl), '.mat'], 'F', 'LABELL')


%% F_imf
% 
% for i = 1:size(SLimf, 2)
%     for j = 1:size(SLimf, 3)
%         shannon_entropy(i, j)   = wentropy(SLimf(:, i, j), 'shannon');
%         mean_frequency(i, j)    = meanfreq(SLimf(:, i, j), fs);
%         spectral_entropy(i, j)  = sum(pentropy(SLimf(:, i, j), fs));
%         tmp                     = aryule(SLimf(:, i, j), 3);
%         AR_coeff(i, j, 1:3)     = tmp(2:4);
%     end
% end
% [M, N, O] = size(AR_coeff);
% AR_coeff    = reshape(AR_coeff, M, O*N);
% 
% F_SLimf = [squeeze(mean(SLimf)), ... % mean              
%         squeeze(var(SLimf)), ...  % variance
%         squeeze(mode(SLimf)), ... % mode
%         squeeze(median(SLimf)),...   % median
%         squeeze(kurtosis(SLimf)),... % kurotsis
%         squeeze(skewness(SLimf)),... % skewness 
%         squeeze(std(SLimf)),...    % standard deviation
%         squeeze(std(SLimf))./squeeze(mean(SLimf)), ... % coefficient of variation
%         squeeze(quantile(SLimf, .25)),... % first quartile
%         squeeze(quantile(SLimf, .75)),... % third quartile
%         squeeze(iqr(SLimf)),...    % interquartile range
%         squeeze(rms(SLimf)),...    % root mean square
%         squeeze(rms(SLimf).^2),... % average power
%         squeeze(sum(SLimf.^2)),... % energy
%         shannon_entropy, ...    % shannon entropy
%         mean_frequency,...      % mean frequency
%         spectral_entropy,...    % spectral entropy
%         AR_coeff];              % AR coefficients           
% 
% 
% clear shannon_entropy mean_frequency spectral_entropy tmp AR_coeff
% 
% for i = 1:size(SRimf, 2)
%     for j = 1:size(SRimf, 3)
%         shannon_entropy(i, j)   = wentropy(SRimf(:, i, j), 'shannon');
%         mean_frequency(i, j)    = meanfreq(SRimf(:, i, j), fs);
%         spectral_entropy(i, j)  = sum(pentropy(SRimf(:, i, j), fs));
%         tmp                     = aryule(SRimf(:, i, j), 3);
%         AR_coeff(i, j, 1:3)     = tmp(2:4);
%     end
% end
% [M, N, O] = size(AR_coeff);
% AR_coeff    = reshape(AR_coeff, M, O*N);
% 
% F_SRimf = [squeeze(mean(SRimf)), ... % mean              
%         squeeze(var(SRimf)), ...  % variance
%         squeeze(mode(SRimf)), ... % mode
%         squeeze(median(SRimf)),...   % median
%         squeeze(kurtosis(SRimf)),... % kurotsis
%         squeeze(skewness(SRimf)),... % skewness 
%         squeeze(std(SRimf)),...    % standard deviation
%         squeeze(std(SRimf))./squeeze(mean(SRimf)), ... % coefficient of variation
%         squeeze(quantile(SRimf, .25)),... % first quartile
%         squeeze(quantile(SRimf, .75)),... % third quartile
%         squeeze(iqr(SRimf)),...    % interquartile range
%         squeeze(rms(SRimf)),...    % root mean square
%         squeeze(rms(SRimf).^2),... % average power
%         squeeze(sum(SRimf.^2)),... % energy
%         shannon_entropy, ...    % shannon entropy
%         mean_frequency,...      % mean frequency
%         spectral_entropy,...    % spectral entropy
%         AR_coeff];              % AR coefficients   
% 
%     
% clear shannon_entropy mean_frequency spectral_entropy tmp AR_coeff
% 
% for i = 1:size(TRimf, 2)
%     for j = 1:size(TRimf, 3)
%         shannon_entropy(i, j)   = wentropy(TRimf(:, i, j), 'shannon');
%         mean_frequency(i, j)    = meanfreq(TRimf(:, i, j), fs);
%         spectral_entropy(i, j)  = sum(pentropy(TRimf(:, i, j), fs));
%         tmp                     = aryule(TRimf(:, i, j), 3);
%         AR_coeff(i, j, 1:3)     = tmp(2:4);
%     end
% end
% [M, N, O] = size(AR_coeff);
% AR_coeff    = reshape(AR_coeff, M, O*N);
% 
% F_TRimf = [squeeze(mean(TRimf)), ... % mean              
%         squeeze(var(TRimf)), ...  % variance
%         squeeze(mode(TRimf)), ... % mode
%         squeeze(median(TRimf)),...   % median
%         squeeze(kurtosis(TRimf)),... % kurotsis
%         squeeze(skewness(TRimf)),... % skewness 
%         squeeze(std(TRimf)),...    % standard deviation
%         squeeze(std(TRimf))./squeeze(mean(TRimf)), ... % coefficient of variation
%         squeeze(quantile(TRimf, .25)),... % first quartile
%         squeeze(quantile(TRimf, .75)),... % third quartile
%         squeeze(iqr(TRimf)),...    % interquartile range
%         squeeze(rms(TRimf)),...    % root mean square
%         squeeze(rms(TRimf).^2),... % average power
%         squeeze(sum(TRimf.^2)),... % energy
%         shannon_entropy, ...    % shannon entropy
%         mean_frequency,...      % mean frequency
%         spectral_entropy,...    % spectral entropy
%         AR_coeff];              % AR coefficients   
% 
% clear shannon_entropy mean_frequency spectral_entropy tmp AR_coeff
% 
% for i = 1:size(TLimf, 2)
%     for j = 1:size(TLimf, 3)
%         shannon_entropy(i, j)   = wentropy(TLimf(:, i, j), 'shannon');
%         mean_frequency(i, j)    = meanfreq(TLimf(:, i, j), fs);
%         spectral_entropy(i, j)  = sum(pentropy(TLimf(:, i, j), fs));
%         tmp                     = aryule(TLimf(:, i, j), 3);
%         AR_coeff(i, j, 1:3)     = tmp(2:4);
%     end
% end
% [M, N, O] = size(AR_coeff);
% AR_coeff  = reshape(AR_coeff, M, O*N);
% 
% F_TLimf = [squeeze(mean(TLimf)), ... % mean              
%         squeeze(var(TLimf)), ...  % variance
%         squeeze(mode(TLimf)), ... % mode
%         squeeze(median(TLimf)),...   % median
%         squeeze(kurtosis(TLimf)),... % kurotsis
%         squeeze(skewness(TLimf)),... % skewness
%         squeeze(std(TLimf)),...    % standard deviation
%         squeeze(std(TLimf))./squeeze(mean(TLimf)), ... % coefficient of variation
%         squeeze(quantile(TLimf, .25)),... % first quartile
%         squeeze(quantile(TLimf, .75)),... % third quartile
%         squeeze(iqr(TLimf)),...    % interquartile range
%         squeeze(rms(TLimf)),...    % root mean square
%         squeeze(rms(TLimf).^2),... % average power
%         squeeze(sum(TLimf.^2)),... % energy
%         shannon_entropy, ...    % shannon entropy
%         mean_frequency,...      % mean frequency
%         spectral_entropy,...    % spectral entropy
%         AR_coeff];              % AR coefficients   
% 
% F_SLimf(:, end+1) =   zeros(size(F_SLimf, 1), 1);
% F_SRimf(:, end+1) =   ones(size(F_SRimf, 1), 1);
% F_TLimf(:, end+1) = 2*ones(size(F_TLimf, 1), 1);
% F_TRimf(:, end+1) = 3*ones(size(F_TRimf, 1), 1);
% 
% Fimf = [F_SLimf; F_SRimf; F_TLimf; F_TRimf];
% 
% kk = 0;
% for i = 0:IMF_n
%     for j = 1:numel(LABELL)
%         kk = kk + 1;
%         if i == 0
%             LABELL2{kk} = LABELL{j};
%         else
%             LABELL2{kk} = [LABELL{j}, '-IMF', num2str(i)];
%         end
%     end
% end
% 
% LABELL = LABELL2;
% 
% save([pwd, '/../Fimf_', num2str(tint), '_', num2str(wl), '.mat'], 'Fimf', 'LABELL')
% 

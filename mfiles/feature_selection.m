
clc
clear

%%

total    = 1;
class    = 2;

MatName  = 'F_30_2.mat';
k        = 10;
TreeCoef = 50;

thr      = 0.315;

%%


load([pwd, '/../', MatName]);


%%

if class == 2

% SL-TL

    label0 = find(F(:, end) == 0);
    label1 = find(F(:, end) == 2);
    Fnew   = F([label0; label1], :);
    Fnew(Fnew(:, end) == 0, end) = 0;
    Fnew(Fnew(:, end) == 2, end) = 1;

    
% TL-TR

%     label2 = find(F(:, end) == 2);
%     label3 = find(F(:, end) == 3);
%     Fnew   = F([label2; label3], :);
%     Fnew(Fnew(:, end) == 2, end) = 0;
%     Fnew(Fnew(:, end) == 3, end) = 1;




    F = Fnew;
    
   
end

%%


rng (1) % For reproducibility

tmp = strfind(MatName, '_');
if total == 1
    if tmp(1) == 2
        FF = F;
    else
        FF = Fimf;
    end
else
    
    load([pwd, '/../fidd']);
    
    if tmp(1) == 2
        FF = F(:, [fidd, size(F, 2)]); 
    else
        FF = Fimf(:, [fidd, size(Fimf, 2)]); 
    end 
end

indx       = randperm(size(FF,1));
X          = FF(indx,1:end-1); %X=a
Y          = FF(indx,end);     %y=t
n          = size(FF,1);  
TTT        = 0;
paroptions = statset('UseParallel', true, 'display', 'final');



%% feature selection



X_train = X(1:round(0.9*length(indx)), :);
Y_train = Y(1:round(0.9*length(indx)), :);

X_test  = X(round(0.9*length(indx))+1:end, :);
Y_test  = Y(round(0.9*length(indx))+1:end, :);

c       = cvpartition(Y_train, 'KFold', k);

classf  = @(train_data, train_labels, test_data, test_labels)...
    sum(predict(fitcsvm(train_data, train_labels), test_data) ~= test_labels);

[fs, history] = sequentialfs(classf, X_train, Y_train, 'cv', c, 'options', paroptions, 'nfeatures', 5);



% ++++++ 30-3 ++++++

% (((SVM)))
% 4 class
% 10   -    [52 128 241 261 341 342 362 370 398 405]; 
% 5    -    [313 341 342 367, 380]
 
% 2 class
% (SL-TL) 
% 5    -    [99 313 315 339 361]
% 10   -    [1 99 116 122 313 315 318 339 353 361]

% (TL-TR) 
% 10   -    [1 2 3 4 141 193 320 345 380 405]
% 5    -    [141 320 345 380 405]



% (((KNN)))
% 4 class
% 10   -    []
% 5    -    []
 
% 2 class
% (SL-TL) 
% 5    -    []
% 10   -    []

% (TL-TR) 
% 10   -    []
% 5    -    []


% ++++++ 30-2 ++++++

% (((SVM)))
% 4 class
% 10   -    [98 241 318 319 341 362 367 370 405 420 ]
% 5    -    [241 318 341 362 405]
 
% 2 class
% (SL-TL) 
% 5    -    []
% 10   -    []

% (TL-TR) 
% 10   -    [1 2 3 4 123 131 319 358 363 383]
% 5    -    [123 131 319 358 383]



% (((KNN)))
% 4 class
% 10   -    []
% 5    -    []
 
% 2 class
% (SL-TL) 
% 5    -    []
% 10   -    []

% (TL-TR) 
% 10   -    []
% 5    -    []



clc
clear

%%

total    = 1;
class    = 4;

MatName  = 'F_30_3.mat';
k        = 10;
TreeCoef = 50;

thr      = 0.315;

%%


load([pwd, '/../', MatName]);


%%

if class == 2
    
    %label0 = find(F(:, end) == 0);
    label1 = find(F(:, end) == 1);
    %label2 = find(F(:, end) == 2);
    label3 = find(F(:, end) == 3);

    Fnew = F([label1; label3], :);
    Fnew(Fnew(:, end) == 1, end) = 0;
    %Fnew(Fnew(:, end) == 2, end) = 0;
    Fnew(Fnew(:, end) == 3, end) = 1;

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


% fs   = [27, 35, 145, 198, 246];
% FF   = F(:, [fs, size(F, 2)]);

fs   = [52 128 241 261 341 342 362 370 398 405]; % 30-3, 10 ta
FF   = F(:, [fs, size(F, 2)]);

indx = randperm(size(FF,1));
X    = FF(indx,1:end-1); %X=a
Y    = FF(indx,end);     %y=t
n    = size(FF,1);  
TTT  = 0;
c    = cvpartition(Y, 'KFold', k);
paroptions = statset('UseParallel', true, 'display', 'final');



%% feature selection


% 
% X_train = X(1:round(0.9*length(indx)), :);
% Y_train = Y(1:round(0.9*length(indx)), :);
% 
% X_test  = X(round(0.9*length(indx))+1:end, :);
% Y_test  = Y(round(0.9*length(indx))+1:end, :);
% 
% c       = cvpartition(Y_train, 'KFold', k);
% 
% classf  = @(train_data, train_labels, test_data, test_labels)...
%     sum(predict(fitcecoc(train_data, train_labels), test_data) ~= test_labels);
% 
% % classf  = @(train_data, train_labels, test_data, test_labels)...
% %     sum(str2double(cell2mat(predict(TreeBagger(50, train_data, train_labels,'OOBPrediction','On', ...
% %                           'Method','classification','OOBPredictorImportance','on', ...
% %                            'SampleWithReplacement','On','InBagFraction',1, ...
% %                            'MinLeafsize',1, 'option', paroptions ...
% %                            ), ...
% %                 test_data))) ~= test_labels);
% 
% [fs, history] = sequentialfs(classf, X_train, Y_train, 'cv', c, 'options', paroptions, 'nfeatures', 10);



%%

for i = 1:k
    
    indxTrn   = c.training(i);    % Training set indices
    indxTest  = c.test(i);        % Test set indices
    

    Mdl = fitcecoc(X(indxTrn,:), Y(indxTrn),'OptimizeHyperparameters','auto',...
          'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
          'expected-improvement-plus','ShowPlots', false)); % Bayes' Optimizatio

%     Mdl       = TreeBagger(TreeCoef, X(indxTrn,:),Y(indxTrn),'OOBPrediction','On','Method','classification','OOBPredictorImportance','on','SampleWithReplacement','On','InBagFraction',1,'MinLeafsize',1, 'option', paroptions);   % 'option', paroptions
    
    label     = predict(Mdl,X(indxTest,:));

    %   acc(i) = (sum(label == Y(indxTest) )/numel(Y(indxTest)));
    
    Tcm{i} = confusionmat(label,double(Y(indxTest)));
    TTT    = TTT + Tcm{i};
%     imp{i} = Mdl.OOBPermutedPredictorDeltaError;
    
    i
end


%%

if class == 2
    
    adad = [length(find(Y==0)), length(find(Y==1)) ;
            length(find(Y==0)), length(find(Y==1))];
    
else
    
    adad = [length(find(Y==0)), length(find(Y==1)), length(find(Y==2)), length(find(Y==3));
            length(find(Y==0)), length(find(Y==1)), length(find(Y==2)), length(find(Y==3));
            length(find(Y==0)), length(find(Y==1)), length(find(Y==2)), length(find(Y==3));
            length(find(Y==0)), length(find(Y==1)), length(find(Y==2)), length(find(Y==3))];

end

TTTper = 100*TTT./adad;

accuracy = 100*sum(diag(TTT))/numel(Y)

if class == 2 
    
    precision   = 100*(TTT(1,1)/sum(TTT(:,1)));
    
    sensitivity = 100*(TTT(1,1)/sum(TTT(1,:)));
    
    specificty  = 100*(TTT(2,2)/sum(TTT(2,:)));

else
    
    precisionSL   = 100*(TTT(1,1)/sum(TTT(:,1)));
    precisionSR   = 100*(TTT(2,2)/sum(TTT(:,2)));
    precisionTL   = 100*(TTT(3,3)/sum(TTT(:,3)));
    precisionTR   = 100*(TTT(4,4)/sum(TTT(:,4)));
    
    sensitivitySL = 100*(TTT(1,1)/sum(TTT(1,:)));
    sensitivitySR = 100*(TTT(2,2)/sum(TTT(2,:)));
    sensitivityTL = 100*(TTT(3,3)/sum(TTT(3,:)));
    sensitivityTR = 100*(TTT(4,4)/sum(TTT(4,:)));
    
    specifictySL  = 100*(sum(sum(TTT(2:end , 2:end))) / sum(sum(TTT(2:end, :))));
    specifictySR  = 100*(sum(sum(TTT([1,3,4],[1,3,4])))/sum(sum(TTT([1,3,4], :))));
    specifictyTL  = 100*(sum(sum(TTT([1,2,4],[1,2,4])))/sum(sum(TTT([1,2,4], :))));
    specifictyTR  = 100*(sum(sum(TTT([1,2,3],[1,2,3])))/sum(sum(TTT([1,2,3], :)))); 
    
end

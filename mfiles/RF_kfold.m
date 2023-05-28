
clc
clear

%%

total    = 0;
class    = 2;

MatName  = 'F_30_2.mat';
k        = 10;
TreeCoef = 50;

thr      = 0.25;

%%


load([pwd, '/../', MatName]);


%%

if class == 2
    
    label0 = find(F(:, end) == 0);
%     label1 = find(F(:, end) == 1);
    label2 = find(F(:, end) == 2);
%     label3 = find(F(:, end) == 3);

    Fnew = F([label0; label2], :);
    Fnew(Fnew(:, end) == 0, end) = 0;
    %Fnew(Fnew(:, end) == 2, end) = 0;
    Fnew(Fnew(:, end) == 2, end) = 1;

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
    
    l1 = find(F(:, end) == 0);
    l2 = find(F(:, end) == 1);
    for i = 1:numel(fidd)
        x  = F(l1, fidd(i));
        y  = F(l2, fidd(i));

        d(i)  = mean(y)-mean(x);


        [hh, pp]  = ttest2(x, y);
        h(i) = hh;
        p(i) = pp;
        fprintf([selected_feature{1, i}, '   : d = %1.3d,    p-val = %1.3d\n'], d(i), p(i));
        
    end
    
    
end


% fs = [342, 341, 333, 314, 384];
% FF = F(:, [fs, size(F, 2)]);








indx = randperm(size(FF,1));
X    = FF(indx,1:end-1); %X=a
Y    = FF(indx,end);     %y=t
n    = size(FF,1);  
TTT  = 0;
c    = cvpartition(Y, 'KFold', k);
paroptions = statset('UseParallel', true, 'display', 'iter');


%%

classs = unique(FF(:, end));

for i = 1:numel(classs)
    sample(i) = numel(find(FF(:, end) == classs(i)));
    gp{i}     = FF(FF(:, end) == classs(i), 1:end-1);
    gp{i}     = gp{i}(randperm(size(gp{i},1)), :);
end

sample = min(sample);

for i = 1:numel(classs)
    gp{i} = gp{i}(1:sample, :);
end


%%

for i = 1:k
    
    indxTrn   = c.training(i);    % Training set indices
    indxTest  = c.test(i);        % Test set indices
    Mdl       = TreeBagger(TreeCoef, X(indxTrn,:),Y(indxTrn),'OOBPrediction','On','Method','classification','OOBPredictorImportance','on','SampleWithReplacement','On','InBagFraction',1,'MinLeafsize',1, 'option', paroptions);   % 'option', paroptions
    label     = predict(Mdl,X(indxTest,:));

    %   acc(i) = (sum(label == Y(indxTest) )/numel(Y(indxTest)));
    
    Tcm{i} = confusionmat(str2num(cell2mat(label)),double(Y(indxTest)));
    TTT    = TTT + Tcm{i};
    imp{i} = Mdl.OOBPermutedPredictorDeltaError;
    
    i
end

% figure
view(Mdl.Trees{50},'Mode','graph')
% F1s(4) = ans;

% meanacc = 100*mean(acc(:));
% meanaccClass= 100*mean(accClass,2);
% imp = oobPermutedPredictorImportance(Mdl);
% imp = OOBPredictorImportance(Mdl);
% imp = Mdl.OOBPermutedPredictorDeltaError;
% idd = impp>0.5;
% fidd = find(idd);

%%

importance = zeros(size(imp{1}));
for i = 1:numel(imp)
    importance = importance + imp{i};
end

importance = importance/numel(imp);
idd  = importance > thr;
fidd = find(idd);

if total == 1
    selected_feature_importance = importance(fidd);
    selected_feature            = LABELL(fidd);

    [selected_feature_importance, I]       = sort(selected_feature_importance, 'descend');
    selected_feature                       = selected_feature(I);
    for i = 1:numel(I)
        selected_feature{2, i} = selected_feature_importance(i);
    end
    
    fidd = fidd(I);

    save([pwd, '/../fidd'], 'fidd', 'selected_feature');
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

TTTper   = 100*TTT./adad;
accuracy = 100*sum(diag(TTT))/numel(Y)

if class == 2 
    
    precision     = 100*(TTT(1,1)/sum(TTT(:,1)));
    sensitivity   = 100*(TTT(1,1)/sum(TTT(1,:)));
    specificty    = 100*(TTT(2,2)/sum(TTT(2,:)));

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



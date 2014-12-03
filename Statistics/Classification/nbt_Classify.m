function [s,ModelVars,Bioms]=nbt_Classify(BCell,Outcome,s,Type, ChannelsOrRegionsToUse)
NCrossVals=20;

%% create DataMatrix from BCell:

DataMatrix = extract_BCell(BCell);

% n_group1=size(BCell{1},2); % no of subjects in the first group
% n_group2=size(BCell{2},2); % no of subjects in the second group

% Outcome = [zeros(1,n_group1) ones(1,n_group2)]'; % n x 1

%% Set up parameters && remove NaN biomarkers
% It's possible to treat the NaN's in another way; interpolate them with EM
% algorithms, set them to 0, etc. Now, we remove ...
% the entire biomarker if it contains any missing data.
%[~,Y]=find(isnan(DataMatrix)==1);
%Y=unique(Y);
%DataMatrix=DataMatrix(:,setdiff(1:size(DataMatrix,2),Y));
%DataMatrix(1:92,:) = DataMatrix(1:92,:)-DataMatrix(93:end,:);
%DataMatrix(93:end,:) = -1.*DataMatrix(1:92,:);
%DataMatrix = sign(DataMatrix).*log10(abs(DataMatrix));


% [c, Sample, l] = princomp(DataMatrix);
% tmp=cumsum(l);
% nr = find(tmp/tmp(end)>0.999,1)
% DataMatrix = Sample(:,1:nr);
% ChannelsToUse =[];

Outcome = Outcome-1;
Outcome = Outcome.';



% TargetNum = 2*length(find(Outcome == 0 ));
% [Outcome, sortIndex]= sort(Outcome);
% DataMatrix = DataMatrix(sortIndex,:);
% DataMatrix = DataMatrix(1:TargetNum,:);
% Outcome = Outcome(1:TargetNum);

%save DataMatrix DataMatrix %sorry Sonja :(
%save Outcome Outcome  %also not saving s further down

% Outcome = Outcome(randperm(length(Outcome)));
 DataMatrix = randn(size(DataMatrix,1),size(DataMatrix,2));
% DataMatrix = DataMatrix(randperm(size(DataMatrix,1)),:);
% warning('random outcome')

Bioms=NaN(size(DataMatrix,2),NCrossVals);
ModelVars=cell(NCrossVals,1);
%% Divide into Train and test set
switch lower(Type)
    case 'crossvalidate'
        % Type CrossValidate
        disp('Cross validation needs work')
        %   DataMatrix = abs(DataMatrix);
        %   DataMatrix = zscore(DataMatrix);
        %       if ~isstruct(ChannelsOrRegionsToUse) && length(ChannelsOrRegionsToUse)>1 % using channels, not regions
        %           [DataMatrix, BiomsToUse] = nbt_RemoveFeatures( DataMatrix,Outcome,'all',ChannelsOrRegionsToUse, size(BCell{1},3));
        %       end
        
        % For this type we randomly
        TestLimit = floor(size(DataMatrix,1)*1/3); %a potential parameter!
        tic
        
        %     BiomsToUse{1,1} = [1:1:length(ChannelsToUse)]'; % use all channels
        
        for i=1:NCrossVals % also potential parametere!
            disp(i)
            try
                
                [ TrainMatrix,  TestMatrix, TrainOutcome, TestOutcome] = ...
                    nbt_RandomSubsampler( DataMatrix,Outcome,TestLimit,'stratified');
                %We use a stratified sample to preserve the class balance. 
                
                %% we first remove features
                [TrainMatrix, BiomsToUse] = nbt_RemoveFeatures( TrainMatrix,TrainOutcome,'ttest2-MCP',ChannelsOrRegionsToUse, size(BCell{1},3));
                
                if(size(BiomsToUse,2) ==1 && size(BCell{1},3) ~= 1)
                    TestMatrix = TestMatrix(:,BiomsToUse{1,1});
                else
                    NewTestMatrix = nan(size(TestMatrix,1),size(BCell{1},3));
                    for ii=1:size(BCell{1},3)
                        NewTestMatrix(:,ii) = nanmedian(TestMatrix(:,BiomsToUse{1,ii}),2);
                    end
                    TestMatrix = NewTestMatrix;
                    TestMatrix=TestMatrix(:,~isnan(TestMatrix(1,:)));
                    clear NewTestMatrix;
                end
                
                %let's also balance the classes
                             TargetNum = 2*length(find(TrainOutcome == 0 ));
                             [TrainOutcome, sortIndex]= sort(TrainOutcome);
                             TrainMatrix = TrainMatrix(sortIndex,:);
                             TrainMatrix = TrainMatrix(1:TargetNum,:);
                             TrainOutcome = TrainOutcome(1:TargetNum);
                
                %             TargetNum = 2*length(find(TestOutcome == 0 ));
                %             TestMatrix = TestMatrix(1:TargetNum,:);
                %             TestOutcome = TestOutcome(1:TargetNum);
                
                %              s.statfunc = 'elasticlogit';
                %              [s] = nbt_TrainClassifier(TrainMatrix,TrainOutcome, s);
                %              if(isempty(find(s.ModelVar(2:end))))
                %                  [s] = nbt_TrainClassifier(TrainMatrix,TrainOutcome, s);
                %              end
                %              TrainMatrix = TrainMatrix(:,find(s.ModelVar(2:end)));
                %              TestMatrix = TestMatrix(:,find(s.ModelVar(2:end)));
                %              s.statfunc = 'lssvm';
                [s] = nbt_TrainClassifier(TrainMatrix,TrainOutcome, s);
                [pp, s ] = nbt_UseClassifier(TestMatrix, s);
                [FPt, TPt, FNt, TNt, SEt, SPt, PPt, NNt, LPt, LNt, MMt, AUCt,H2] = ...
                    nbt_evalOutcome(pp, TestOutcome);
            catch
                FPt = nan;
                TPt = nan;
                FNt = nan;
                TNt = nan;
                SEt = nan;
                SPt = nan;
                PPt = nan;
                NNt = nan;
                LPt = nan;
                LNt = nan;
                MMt = nan;
                AUCt = nan;
                H2 = nan;
                s.ModelVar = [];
            end
            
            
            % training and testing on the same data
            %             [pp, s ] = nbt_UseClassifier(TrainMatrix, s);
            %             [FPt, TPt, FNt, TNt, SEt, SPt, PPt, NNt, LPt, LNt, MMt, AUCt,H2] = ...
            %                 nbt_evalOutcome(pp, TrainOutcome);
            
            ModelVars{i}=s.ModelVar;
            if(iscell(FPt))
                for GrpID = 1:size(FPt,1)
                    FP{GrpID,1}(i) = FPt{GrpID,1};
                    TP{GrpID,1}(i) = TPt{GrpID,1};
                    FN{GrpID,1}(i) = FNt{GrpID,1};
                    TN{GrpID,1}(i) = TNt{GrpID,1};
                    SE{GrpID,1}(i) = SEt{GrpID,1};
                    SP{GrpID,1}(i) = SPt{GrpID,1};
                    PP{GrpID,1}(i) = PPt{GrpID,1};
                    NN{GrpID,1}(i) = NNt{GrpID,1};
                    LP{GrpID,1}(i) = LPt{GrpID,1};
                    LN{GrpID,1}(i) = LNt{GrpID,1};
                    MM{GrpID,1}(i) = MMt{GrpID,1};
                end
            else %FTt is not a cell
                FP(i) =  FPt;
                TP(i) =  TPt;
                FN(i) =  FNt;
                TN(i) =  TNt;
                SE(i) =  SEt;
                SP(i) =  SPt;
                PP(i) =  PPt;
                NN(i) =  NNt;
                LP(i) =  LPt;
                LN(i) =  LNt;
                MM(i) =  MMt;
                AUC(i) = AUCt;
                H_measure(i)=H2;
                accuracy(i) = (TP(i)+TN(i))./(TN(i)+TP(i)+FN(i)+FP(i))*100;
            end
            
        end
        disp('CrossValidate:done')
        toc
    case 'validate'
        tic
        % Type Validate
        TestLimit = floor(size(DataMatrix,1)*1/3); %a potential parameter!
        [ TrainMatrix,  TestMatrix, TrainOutcome, TestOutcome] = ...
            nbt_RandomSubsampler( DataMatrix,Outcome,TestLimit,'simple');
        %We use a stratified sample to preserve the class balance.
        %% Feature selction - we first prune the biomarkers given to the classsification algorithm
        % [TrainMatrix, BiomsToUse] = nbt_RemoveFeatures( TrainMatrix,TrainOutcome,'ttest2',ChannelsToUse, size(BCell{1},3));
        
        
        % call nbt_TrainClassifier
        [s] = nbt_TrainClassifier(TrainMatrix,TrainOutcome, s);
        
        if(size(BiomsToUse,2) ==1 && size(BCell{1},3) ~= 1)
            TestMatrix = TestMatrix(:,BiomsToUse{1,1});
        else
            NewTestMatrix = nan(size(TestMatrix,1),size(BCell{1},3));
            for ii=1:size(BCell{1},3)
                NewTestMatrix(:,ii) = nanmedian(TestMatrix(:,BiomsToUse{1,ii}),2);
            end
            TestMatrix = NewTestMatrix;
            clear NewTestMatrix;
        end
        % call nbt_TestClassifier
        [pp, s] = nbt_UseClassifier(TestMatrix, s);
        %eval outcome
        [FP, TP, FN, TN, SE, SP, PP, NN, LP, LN, MM, AUC, H_measure] = nbt_evalOutcome(pp, TestOutcome);
        toc
    case 'train'
        % Type Train
        s = nbt_TrainClassifier(DataMatrix, Outcome, s);
    case 'use'
        % Type Use
        [pp, s] = nbt_UseClassifier(DataMatrix, s);
        %eval outcome
        [FP, TP, FN, TN, SE, SP, PP, NN, LP, LN, MM] = nbt_evalOutcome(pp, TestOutcome);
    otherwise
        error('The specified Type is not known - Please use either CrossValidate, Validate, Train, or Use')
end

%update s
s.FalsePositive =  FP;
s.TruePositive =  TP;
s.FalseNegative =  FN;
s.TrueNegative =  TN;
s.Sensitivity =  SE;
s.Specificity =  SP;
s.PositivePredictiveValue  =  PP;
s.NegativePredictiveValue =  NN;
s.LikelihoodRatioPos =  LP;
s.LikelihoodRatioNeg  =  LN;
s.MatthewCorr =  MM;
s.AUC=AUC;
s.H_measure=H_measure;


%save s s


%                 s.BaselineSE=BaselineSE;
%                 s.BaselineSP=BaselineSP;
%                 s.BaselineAUC=BaselineAUC;
%                 s.ESE=ESE;
%                 s.ESP=ESP;
%                 s.EAUC=EAUC;
%                 s.HAE=HAE;
%                 s.HE=HE;

%% plotting
%first calculate pp for the full matrix
[pp] = nbt_UseClassifier(DataMatrix, s);
% make pp values for each group
pp1 = pp(Outcome == 0);
pp2 = pp(Outcome == 1);
% DataMatrix = zscore(DataMatrix)';
DataMatrix = DataMatrix';
DataReal1 = DataMatrix(:, Outcome == 0);
DataReal2 = DataMatrix(:, Outcome == 1);
DataPred1 = DataMatrix(:, pp < 0.5);
DataPred2 = DataMatrix(:, pp >= 0.5);

figure('Name','Classification performance')

% plot dot plot of pp values for the full Data matrix
nbt_DotPlot(subplot(2,4,1),0.1,0.025,0,0,@median,{'Group 1';'Group 2'; 'Probability'},'',pp1',1:length(pp1),1,pp2',1:length(pp2),1);
set(gca,'YLim',[0 1])
% plot dot plot of Data values - real and classified
%first real
nbt_DotPlot(subplot(2,4,2),0.1,0.025,0,0,@median,{'Group 1';'Group 2'; 'Biomarker values'},'',DataReal1,1:size(DataReal1,2),1:size(DataReal1,1),DataReal2,1:size(DataReal2,2),1:size(DataReal2,1));
% then predicted
nbt_DotPlot(subplot(2,4,3),0.1,0.025,0,0,@median,{'Group 1';'Group 2'; 'Predicted groups: Biomarker values'},'',DataPred1,1:size(DataPred1,2),1:size(DataPred1,1),DataPred2,1:size(DataPred2,2),1:size(DataPred2,1));
% Plot ROC
subplot(2,4,4)
[FPR,TPR] = perfcurve(Outcome,pp,1);
plot(FPR,TPR)
xlabel('False positive rate'); ylabel('True positive rate')
title('ROC')

subplot(2,4,5)
boxplot(accuracy/100)
set(gca,'YLim',[0.5 1])
title('Accuracy')

subplot(2,4,6)
boxplot(SE)
set(gca,'YLim',[0.5 1])
title('Sensitivity (SE)')

subplot(2,4,7)
boxplot(SP)
set(gca,'YLim',[0.5 1])
title('Specificity (SP)')

subplot(2,4,8)
boxplot(PP)
set(gca,'YLim',[0.5 1])
title('Precision (PP)')


%% nested function part
    function DataMatrix = extract_BCell(BCell)
        if(isempty(ChannelsOrRegionsToUse))
            ChannelsOrRegionsToUse = 1:size(BCell{1},1);
        end
        
        for ii=1:size(BCell{1},3)
            disp(ii)
            if isstruct(ChannelsOrRegionsToUse)
                ChansOrRegsToUse = [1: size(ChannelsOrRegionsToUse,2)]; % all regions
            else
                ChansOrRegsToUse = ChannelsOrRegionsToUse; % single region
            end
            
            for i=ChansOrRegsToUse;
                if ~exist('becell');
                    becell{1,1}=BCell{1}(i,:,ii);
                    becell{2,1}=BCell{2}(i,:,ii);
                else
                    becell{1,1}(end+1,:,:)=BCell{1}(i,:,ii);
                    becell{2,1}(end+1,:,:)=BCell{2}(i,:,ii);
                end
            end
        end
        for GrpID = 1:size(becell,1)
            GroupSize(GrpID) = size(becell{GrpID,1},2);
        end
        DataIndex = 1;
        for GrpID = 1:size(becell,1)
            if isrow(squeeze(reshape(becell{GrpID,1}(:,:,:),1,size(becell{GrpID,1},2),size(becell{GrpID,1},1)*size(becell{GrpID,1},3))));
                DataMatrix(DataIndex:DataIndex+GroupSize(GrpID)-1,:) = ...
                    becell{GrpID,1}(:,:,:).';
            else
                DataMatrix(DataIndex:DataIndex+GroupSize(GrpID)-1,:) = becell{GrpID,1}(:,:,:).';
            end
            Outcome(DataIndex:DataIndex+GroupSize(GrpID)-1) = GrpID;
            DataIndex = DataIndex+GroupSize(GrpID);
        end
    end
end

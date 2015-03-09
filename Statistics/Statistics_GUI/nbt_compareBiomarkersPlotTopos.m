function nbt_compareBiomarkersPlotTopos(d1,d2,G,B_values1,B_values2,bioms1,bioms2,biomindex,Pvalues,rho,noGroups,splitType,splitValue,regs_or_chans_name)

pos_cursor_unitfig = get(gca,'currentpoint');


clear Pvalues


for n = 1:size(B_values2,1);
    for m = 1: size(B_values1,1);
        [rho(n,m),Pvalues(n,m)] = corr(B_values1(m,:)',B_values2(n,:)','type','Spearman');
    end
end

if strcmp(bioms1{1},'rsq.Answers')
    question = round(pos_cursor_unitfig(1,2));
    Pvalues = Pvalues';
else
    question = round(pos_cursor_unitfig(1,1));
end
global rsq;
if strcmp(bioms1{1},'rsq.Answers')
    figure('Name',['Topoplots: P-values and Rho values for Question  ' num2str(question) '. ''' rsq.Questions{question} ''' and Biomarker ' bioms2{1}],'NumberTitle','off')
else
    figure('Name',['Topoplots: P-values and Rho values for Question  ' num2str(question) '. ' rsq.Questions{question} ' and Biomarker ' bioms1{1}],'NumberTitle','off')
end

set(gcf,'position',[10          80       700      500])
coolWarm = load('nbt_CoolWarm.mat','coolWarm');
coolWarm = coolWarm.coolWarm;
colormap(coolWarm);
subplot(2,2,3)
if strcmp(regs_or_chans_name,'Regions')
    nbt_plot_subregions_hack(log10(Pvalues(question,:)),-2.6,0);
else
    topoplot(log10(double(Pvalues(question,:))), G(1).chansregs.chanloc,'headrad','rim');
end
minPValue = -2.6;% Plot log10(P-Values) to trick colour bar
maxPValue = -0;
cbh = colorbar('westOutside');
caxis([minPValue maxPValue])
set(cbh,'XTick',[-2 -1.3010 -1 0])
set(cbh,'XTicklabel',[0.01 0.05 0.1 1]) %(log scale)
title('P-values')
subplot(2,2,4)
if strcmp(regs_or_chans_name,'Regions')
    mxRho = max(abs(rho(question,:)));
    nbt_plot_subregions_hack((rho(question,:)),-mxRho,mxRho);
else
    topoplot((rho(question,:)), G(1).chansregs.chanloc,'headrad','rim');
end

colorbar('westoutside')
title('Rho')
subplot(2,2,1)
axis off
if strcmp(bioms1{1},'rsq.Answers')
    
    textThis = sprintf(['P-values and Rho values for Question ' num2str(question) '. ' rsq.Questions{question} ' and Biomarker ' regexprep(bioms2{1},'_',' ')]);
    
else
    textThis = sprintf(['P-values and Rho values for Question ' num2str(question) '. ' rsq.Questions{question} ' and Biomarker ' regexprep(bioms1{1},'_',' ')]);
end
nbt_split_title([1 0.5],textThis,200,11);



clear Pvalues

for n = 1:size(B_values2,1);
    for m = 1: size(B_values1,1);
        if splitType == 1
            %split on %
            nns = nnz(isnan(B_values2(n,:)));
            if nns > 0
                quarts = floor((splitValue/100)*(size(B_values2,2)- nns));
                [valk,inds] = sort(B_values2(n,:));
                Group2{n} = inds(1:quarts);
                Group1{n} = inds(1+end-quarts-nns:end-nns);
            else
                quarts = floor((splitValue/100)*size(B_values2,2));
                [indef,inds] = sort(B_values2(n,:));
                Group2{n} = inds(1:quarts);
                Group1{n} = inds(1+end-quarts:end);
            end
            
            if strcmp(bioms1,'rsq.Answers')
                
                B1(n,:) = B_values1(m,Group1{n}); % biomarker for yes group
                B2(n,:) = B_values1(m,Group2{n}); % biomarker for no group
                [indef,Pvalues(n,m),indef] = ttest2(B1(n,:),B2(n,:));
                rho(n,m) = mean(B1(n,:))-mean(B2(n,:));
                MeanB1(n,m) = mean(B_values1(m,Group1{n})); % biomarker for yes group
                MeanB2(n,m) = mean(B_values1(m,Group2{n})); % biomarker for no group
            else
                B1 = B_values1(m,Group1{n}); % biomarker for yes group
                B2 = B_values1(m,Group2{n}); % biomarker for no group
                [indef,Pvalues(n,m),indef] = ttest2(B1,B2);
                rho(n,m) = mean(B1)-mean(B2);
                MeanB1(n,m) = mean(B_values1(m,Group1{n})); % biomarker for yes group
                MeanB2(n,m) = mean(B_values1(m,Group2{n})); % biomarker for no group
            end
        else
            %split on value
            Group1{n} = find(B_values2(n,:)>splitValue);% yes
            Group2{n} = find(B_values2(n,:)<splitValue);% no
            
            if strcmp(bioms1,'rsq.Answers')
                B1 = B_values1(m,Group1{n}); % biomarker for yes group
                B2 = B_values1(m,Group2{n}); % biomarker for no group
                MeanB1(n,m) = mean(B_values1(m,Group1{n})); % biomarker for yes group
                MeanB2(n,m) = mean(B_values1(m,Group2{n})); % biomarker for no group
                [indef,Pvalues(n,m),indef] = ttest2(B1,B2);
                rho(n,m) = mean(B1)-mean(B2);
            else
                B1 = B_values1(m,Group1{n}); % biomarker for yes group
                B2 = B_values1(m,Group2{n}); % biomarker for no group
                [indef,Pvalues(n,m),indef] = ttest2(B1,B2);
                rho(n,m) = mean(B1)-mean(B2);
                MeanB1(n,m) = mean(B_values1(m,Group1{n})); % biomarker for yes group
                MeanB2(n,m) = mean(B_values1(m,Group2{n})); % biomarker for no group
            end
        end
    end
end


if strcmp(bioms1{1},'rsq.Answers')
    
    p =  Pvalues(:,question);
    diffmean = MeanB1(:,question) - MeanB2(:,question);
    minCL = min(min([MeanB1(:,question) MeanB2(:,question)]));
    maxCL = max(max([MeanB1(:,question) MeanB2(:,question)]));
    MeanB1 = MeanB1(:,question);
    MeanB2 = MeanB2(:,question);
    
else
    p =  Pvalues(question,:);
    diffmean = MeanB1(question,:) - MeanB2(question,:);
    minCL = min(min([MeanB1(question,:) MeanB2(question,:)]));
    maxCL = max(max([MeanB1(question,:) MeanB2(question,:)]));
    MeanB1 = MeanB1(question,:);
    MeanB2 = MeanB2(question,:);
    
end

if strcmp(bioms1{1},'rsq.Answers')
    figure('Name',['Unpaired ttest between Top and bottom group split on ' regexprep(bioms2{biomindex},'_',' ') ' Compared with Question ' num2str(question) '. ' rsq.Questions{question}]);
else
    figure('Name',['Unpaired ttest between Top and bottom group split on Question ' num2str(question) '. ' rsq.Questions{question} ' Compared with biomarker ' regexprep(bioms1{biomindex},'_',' ')]);
end



coolWarm = load('nbt_CoolWarm.mat','coolWarm');
coolWarm = coolWarm.coolWarm;
colormap(coolWarm);
subplot(3,3,4)
if strcmp(regs_or_chans_name,'Regions')
    nbt_plot_subregions_hack(MeanB1,minCL,maxCL);
else
    topoplot( MeanB1, G(1).chansregs.chanloc,'headrad','rim');
end
if splitType == 1
    title(['Mean Group Top ' num2str(splitValue) '%']);
else
    title(['Mean Group above ' num2str(splitValue)]);
end
cbh = colorbar('westOutside');
caxis([minCL maxCL])
subplot(3,3,7)
if strcmp(regs_or_chans_name,'Regions')
    nbt_plot_subregions_hack(MeanB2,minCL,maxCL);
else
    topoplot(MeanB2, G(1).chansregs.chanloc,'headrad','rim');
end

caxis([minCL maxCL])
if splitType == 1
    title(['Mean Group Bottom ' num2str(splitValue) '%']);
else
    title(['Mean Group Below ' num2str(splitValue)]);
end
cbh = colorbar('westOutside');
subplot(3,3,[5 8])
if strcmp(regs_or_chans_name,'Regions')
    nbt_plot_subregions_hack(log10(p),-2.6,-0);
else
    topoplot(log10(p), G(1).chansregs.chanloc,'headrad','rim');
end



minPValue = -2.6;% Plot log10(P-Values) to trick colour bar
maxPValue = -0;
cbh = colorbar('westOutside');
caxis([minPValue maxPValue])
set(cbh,'yTick',[-2 -1.3010 -1 0])
set(cbh,'yTicklabel',[0.01 0.05 0.1 1]) %(log scale)
title('P-values for unpaired ttest')
subplot(3,3,[6 9])
minDCL = -(max(abs(diffmean)));
maxDCL = -minDCL;
if strcmp(regs_or_chans_name,'Regions')
    nbt_plot_subregions_hack(diffmean,minDCL,maxDCL);
else
    topoplot(diffmean, G(1).chansregs.chanloc,'headrad','rim');
end
title('Difference of the means of the 2 groups')
cbh = colorbar('westOutside');
subplot(3,3,1)
axis off


if strcmp(bioms1{1},'rsq.Answers')
    textThis = sprintf(['Unpaired ttest between Top and bottom group split on ' regexprep(bioms2{biomindex},'_',' ') ' Compared with Question ' num2str(question) '. ' rsq.Questions{question}]);
else
    textThis = sprintf(['Unpaired ttest between Top and bottom group split on Question ' num2str(question) '. ' rsq.Questions{question} ' Compared with biomarker ' regexprep(bioms1{biomindex},'_',' ')]);
end

nbt_suptitle(textThis);
%nbt_split_title([1 0.5],textThis,200,11);


end
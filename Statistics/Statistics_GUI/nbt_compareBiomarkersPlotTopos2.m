    function nbt_compareBiomarkersPlotTopos2(d1,d2,G1,G2,B_values1,B_values2,bioms1,bioms2,biomindex,Pvalues,rho)
        
        pos_cursor_unitfig = get(gca,'currentpoint');
        question = round(pos_cursor_unitfig(1,1));
        
        figure('Name',['Topoplots: P-values and Rho values for Question ' num2str(question) ' and Biomarker ' regexprep(bioms2{biomindex},'_',' ')],'NumberTitle','off')
        set(gcf,'position',[10          80       700      500])
        coolWarm = load('nbt_CoolWarm.mat','coolWarm');
        coolWarm = coolWarm.coolWarm;
        colormap(coolWarm);
        subplot(2,2,3)
        topoplot(log10(Pvalues(question,:)), G1.chansregs.chanloc,'headrad','rim');
        minPValue = -2;% Plot log10(P-Values) to trick colour bar
        maxPValue = -0.5;
        cbh = colorbar('westOutside');
        caxis([minPValue maxPValue])
        set(cbh,'XTick',[-2 -1.3010 -1 0])
        set(cbh,'XTicklabel',[0.01 0.05 0.1 1]) %(log scale)
        title('P-values')
        subplot(2,2,4)
        topoplot((rho(question,:)), G1.chansregs.chanloc,'headrad','rim');
        colorbar('westoutside')
        title('Rho')
        subplot(2,2,1)
        axis off
        textThis = sprintf(['P-values and Rho values for Question ' num2str(question) ' and Biomarker ' regexprep(bioms2{biomindex},'_',' ')]);
        nbt_split_title([1 0.5],textThis,200,11);
        
        for i = 1: size(B_values1,1)
            Group1{i} = find(B_values1(i,:)>0);% pos
            Group2{i} = find(B_values1(i,:)<0);% neg
        end
        B1 = B_values2(:,Group1{question},biomindex); % biomarker for yes group
        B2 = B_values2(:,Group2{question},biomindex); % biomarker for no group
        for m = 1:size(B_values2,1)
            [h,p(m),c(m,:)] = ttest2(B1(m,:),B2(m,:));
            diffmean(m) = mean(B1(m,:))-mean(B2(m,:));
        end
        MeanB1 = mean(B1,2);
        MeanB2 = mean(B2,2);
        figure('Name',['Unpaired ttest between yes and no group for ' regexprep(bioms2{biomindex},'_',' ')])
        coolWarm = load('nbt_CoolWarm.mat','coolWarm');
        coolWarm = coolWarm.coolWarm;
        colormap(coolWarm);
        subplot(3,3,4)
        topoplot( MeanB1, G1.chansregs.chanloc,'headrad','rim');
        title('Mean Group "Yes"')
        cbh = colorbar('westOutside');
        subplot(3,3,7)
        topoplot(MeanB2, G1.chansregs.chanloc,'headrad','rim');
        title('Mean Group "No"')
        cbh = colorbar('westOutside');
        subplot(3,3,[5 8])
        topoplot(log10(p), G1.chansregs.chanloc,'headrad','rim');
        minPValue = -2;% Plot log10(P-Values) to trick colour bar
        maxPValue = -0.5;
        cbh = colorbar('westOutside');
        caxis([minPValue maxPValue])
        set(cbh,'yTick',[-2 -1.3010 -1 0])
        set(cbh,'yTicklabel',[0.01 0.05 0.1 1]) %(log scale)
        title('P-values for unpaired ttest')
        subplot(3,3,[6 9])
        topoplot(diffmean, G1.chansregs.chanloc,'headrad','rim');
        title('Difference of the means of the 2 groups')
        cbh = colorbar('westOutside');
        subplot(3,3,1)
        axis off
        textThis = sprintf(['Correlation between two groups for Answers yes and no ']);
        nbt_split_title([1 0.5],textThis,200,11);
        
        
    end
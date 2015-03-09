function nbt_compareBiomarkersPlotChansComp2(d1,d2,B_values1,B_values2,bioms1,bioms2,biomindex)
        pos_cursor_unitfig = get(gca,'currentpoint');
        question = round(pos_cursor_unitfig(1,1));
        chan_or_reg = round(pos_cursor_unitfig(1,2));
        if question>0 && question<=size(B_values1,1) && chan_or_reg>0 && chan_or_reg<=size(B_values2,1)
            B = B_values2(chan_or_reg,:,biomindex);
            B_Ans = B_values1(question,:);
            %-----
            figure('Name',['Least-squares fit for Question ' num2str(question) ' and Channel ' num2str(chan_or_reg) ' of ' regexprep(bioms2{biomindex},'_',' ')],'NumberTitle','off')
            set(gcf,'position',[10          80       700      500])
            plot(B_Ans,B,'.')
            lsline
            xlabel(['Answers to Question ' num2str(question)])
            ylabel([regexprep(bioms2{biomindex},'_',' ') ' for Channel ' num2str(chan_or_reg)])
            textThis = sprintf(['Least-squares fit for Question ' num2str(question) ' and Channel ' num2str(chan_or_reg) ' of ' regexprep(bioms2{biomindex},'_',' ')]);
            nbt_split_title([2 0.9],textThis,200,11);
            %----
            for i = 1: size(B_values1,1)
                Group1{i} = find(B_values1(i,:)>0);% yes
                Group2{i} = find(B_values1(i,:)<0);% no
            end
            B1 = B(Group1{question}); % biomarker for yes group
            B2 = B(Group2{question}); % biomarker for no group
            [h,p,c] = ttest2(B1,B2);
            diffmean = mean(B1)-mean(B2);
            figure('Name',['Mean of the Yes and No group for ' regexprep(bioms2{biomindex},'_',' ') ' and Channel ' num2str(chan_or_reg)],'NumberTitle','off')
            set(gcf,'position',[10          80       700      500])
            bar([mean(B1) mean(B2)])
            ylabel(['Mean ' regexprep(bioms2{biomindex},'_',' ') ' in Channel ' num2str(chan_or_reg)])
            xlim([0 3])
            set(gca,'xtick',[0 1 2 3 ],'Xticklabel',{''; ['"Yes" Gorup (n = ' num2str(length(B1)) ')']; ['"No" Gorup (n = ' num2str(length(B2)) ')'];''})
            textThis = sprintf(['Mean of the Yes and No group for ' regexprep(bioms2{biomindex},'_',' ') ' and Channel ' num2str(chan_or_reg)]);
            nbt_split_title([1.5 0.7],textThis,200,11);
            
        end
    end
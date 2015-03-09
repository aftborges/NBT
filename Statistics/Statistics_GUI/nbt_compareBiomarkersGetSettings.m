function nbt_compareBiomarkersGetSettings(d1,d2,ListBiom1, ListBiom2, ListRegion, ListGroup,ListTest,ListDisplay,ListSplit,ListSplitValue,G)
disp('Computing biomarkers comparison ...')
global rsq %set in nbt_compareBiomarkers
% --- get statistics test (one)
bioms_ind1 = get(ListBiom1,'Value');
bioms_name1 = get(ListBiom1,'String');
bioms_name1 = bioms_name1(bioms_ind1);
% --- get channels or regions (one)
regs_or_chans_index = get(ListRegion,'Value');
regs_or_chans_name = get(ListRegion,'String');
regs_or_chans_name = regs_or_chans_name(regs_or_chans_index);
% --- get biomarkers (one or more)
bioms_ind2 = get(ListBiom2,'Value');
bioms_name2 = get(ListBiom2,'String');
bioms_name2 = bioms_name2(bioms_ind2);
% --- get group (one or more)
group_ind = get(ListGroup,'Value');
group_name = get(ListGroup,'String');
group_name = group_name(group_ind);
% --- get Test
test_ind = get(ListTest,'Value');

% --- get Display
display_ind = get(ListDisplay,'Value');


splitType = get(ListSplit,'Value');
splitValue = str2num(get(ListSplitValue,'String'));

group_diffexist = 0;
%%  If more than one group then has to be paired and make sure that
%  both have the same subjects.
if length(group_ind) == 1
    
    Group = G(group_ind);
         if ~isempty(Group.group_difference)
            group_diffexist = 1;
            path = Group.fileslist.path;
            n_files = length(Group.fileslist);
            [B_values_cell1,Sub,Proj,unit] = nbt_checkif_groupdiff(Group,G,n_files,bioms_name1,path);
            [B_values_cell2,Sub,Proj,unit] = nbt_checkif_groupdiff(Group,G,n_files,bioms_name2,path);
            group_diffexist;
            NCHANNELS = length(Group.chansregs.chanloc);
            [dimens_diff1,biomPerChans1,IndexbiomNotPerChans1] = nbt_dimension_check(B_values_cell1,NCHANNELS);
            
            [dimens_diff2,biomPerChans2,IndexbiomNotPerChans2] = nbt_dimension_check(B_values_cell2,NCHANNELS);
            if ~isempty(biomPerChans1)
                B_values1 = nbt_extractBiomPerChans(biomPerChans1,B_values_cell1);
            else
                dimBio = IndexbiomNotPerChans1{dim};
                B_values1 = nbt_extractBiomNotPerChans(dimBio,B_values_cell1);
            end
            
            if ~isempty(biomPerChans2)
                B_values2 = nbt_extractBiomPerChans(biomPerChans2,B_values_cell2);
            else
                
                dimBio = IndexbiomNotPerChans2{1};
                B_values2 = nbt_extractBiomNotPerChans(dimBio,B_values_cell2);
                
            end
            
            bioms1 = bioms_name1;
            bioms2 = bioms_name2;
            
            
        else
            group_diffexist = 0;
            path = Group.fileslist.path;
            n_files = length(Group.fileslist);
            
            bioms1 = bioms_name1;
            bioms2 = bioms_name2;
            % load biomarkers
            for j = 1:n_files % subject
                for l = 1:length(bioms2) % biomarker
                    namefile = Group.fileslist(j).name;
                    biomarker1 = bioms1{1};
                    biomarker2 = bioms2{l};
                    [B_values1(:,j),Sub,Proj,unit1{j}]=nbt_load_analysis(path,namefile,biomarker1,@nbt_get_biomarker,[],[],[]);
                    [B_values2(:,j,l),Sub,Proj,unit2{j,l}]=nbt_load_analysis(path,namefile,biomarker2,@nbt_get_biomarker,[],[],[]);
                end
            end
        end
        
    elseif length(group_ind) == 2
        Grop1 = G(group_ind(1));
        Grop2 = G(group_ind(2));
        bioms1 = bioms_name1;
            bioms2 = bioms_name2;
            
        if ~isempty(Grop1.group_difference)
            
            path = Grop1.fileslist.path;
            n_files = length(Grop1.fileslist);
            [B_values_cell1,Sub,Proj,unit] = nbt_checkif_groupdiff(Grop1,G,n_files,bioms_name1,path);
            [B_values_cell2,Sub,Proj,unit] = nbt_checkif_groupdiff(Grop1,G,n_files,bioms_name2,path);
            
            NCHANNELS = length(Grop1.chansregs.chanloc);
            [dimens_diff1,biomPerChans1,IndexbiomNotPerChans1] = nbt_dimension_check(B_values_cell1,NCHANNELS);
            
            [dimens_diff2,biomPerChans2,IndexbiomNotPerChans2] = nbt_dimension_check(B_values_cell2,NCHANNELS);
            if ~isempty(biomPerChans1)
                B1_values1 = nbt_extractBiomPerChans(biomPerChans1,B_values_cell1);
            else
                dimBio = IndexbiomNotPerChans1{1};
                B1_values1 = nbt_extractBiomNotPerChans(dimBio,B_values_cell1);
            end
            
            if ~isempty(biomPerChans2)
                B1_values2 = nbt_extractBiomPerChans(biomPerChans2,B_values_cell2);
            else
                
                dimBio = IndexbiomNotPerChans2{1};
                B1_values2 = nbt_extractBiomNotPerChans(dimBio,B_values_cell2);
            end
            
        else
            bioms1 = bioms_name1;
            nameG1 = Grop1.fileslist.group_name;
            path1 = Grop1.fileslist.path;
            n_files1 = length(Grop1.fileslist);
            % load biomarkers
            for j = 1:n_files1 % subject
                for l = 1:length(bioms2) % biomarker
                    namefile = Grop1.fileslist(j).name;
                    biomarker1 = bioms1{1};
                    biomarker2 = bioms2{l};
                    [B1_values1(:,j),Sub,Proj,unit1{j}]=nbt_load_analysis(path1,namefile,biomarker1,@nbt_get_biomarker,[],[],[]);
                    [B1_values2(:,j,l),Sub,Proj,unit2{j,l}]=nbt_load_analysis(path1,namefile,biomarker2,@nbt_get_biomarker,[],[],[]);
                    
                end
            end
            
            
        end
        
        if ~isempty(Grop2.group_difference)
            bioms1 = bioms_name1;
            bioms2 = bioms_name2;
            
            path = Grop2.fileslist.path;
            n_files = length(Grop2.fileslist);
            [B_values_cell1,Sub,Proj,unit] = nbt_checkif_groupdiff(Grop2,G,n_files,bioms_name1,path);
            [B_values_cell2,Sub,Proj,unit] = nbt_checkif_groupdiff(Grop2,G,n_files,bioms_name2,path);
            
            NCHANNELS = length(Grop2.chansregs.chanloc);
            [dimens_diff1,biomPerChans1,IndexbiomNotPerChans1] = nbt_dimension_check(B_values_cell1,NCHANNELS);
            
            [dimens_diff2,biomPerChans2,IndexbiomNotPerChans2] = nbt_dimension_check(B_values_cell2,NCHANNELS);
            if ~isempty(biomPerChans1)
                B2_values1 = nbt_extractBiomPerChans(biomPerChans1,B_values_cell1);
            else
                dimBio = IndexbiomNotPerChans1{1};
                B2_values1 = nbt_extractBiomNotPerChans(dimBio,B_values_cell1);
            end
            
            if ~isempty(biomPerChans2)
                B2_values2 = nbt_extractBiomPerChans(biomPerChans2,B_values_cell2);
            else
                
                dimBio = IndexbiomNotPerChans2{1};
                B2_values2 = nbt_extractBiomNotPerChans(dimBio,B_values_cell2);
            end
            
        else
            bioms2 = bioms_name2;
            nameG2 = Grop2.fileslist.group_name;
            path2 = Grop2.fileslist.path;
            n_files2 = length(Grop2.fileslist);
            for j = 1:n_files2 % subject
                for l = 1:length(bioms2) % biomarker
                    namefile = Grop2.fileslist(j).name;
                    biomarker1 = bioms1{1};
                    biomarker2 = bioms2{l};
                    [B2_values1(:,j),Sub,Proj,unit1{j}]=nbt_load_analysis(path2,namefile,biomarker1,@nbt_get_biomarker,[],[],[]);
                    [B2_values2(:,j,l),Sub,Proj,unit2{j,l}]=nbt_load_analysis(path2,namefile,biomarker2,@nbt_get_biomarker,[],[],[]);
                end
            end
            
            
        end
        B_values1 = B1_values1-B2_values1;
        B_values2 = B1_values2-B2_values2;
        
    end
    
    
    
    if strcmp(regs_or_chans_name,'Regions')
        if strcmp(bioms_name1,'rsq.Answers') == 0
            regions = G(group_ind(1)).chansregs.listregdata;
            for j = 1:size(B_values1,2) % subject
                
                B1 = B_values1(:,j);
                B_gebruik1(:,j) = nbt_compare_getRegions(B1,regions);
            end
            clear B_values1;
            B_values1 = B_gebruik1;
        end
        if strcmp(bioms_name2,'rsq.Answers') == 0
            regions = G(group_ind(1)).chansregs.listregdata;
            for j = 1:size(B_values2,2) % subject
                B2 = B_values2(:,j);
                B_gebruik2(:,j) = nbt_compare_getRegions(B2,regions);
            end
            clear B_values2;
            B_values2 = B_gebruik2;
        end
    end
    
    %%test_ind == 1 is spearman correlation
    %%Otherwise its a ttest and so data needs to be split
    
    if (test_ind == 1)
        for i = 1:size(B_values1,3)
            for n = 1:size(B_values2,1);
                for m = 1: size(B_values1,1);
                    [rho(n,m),Pvalues(n,m)] = corr(B_values1(m,:)',B_values2(n,:)','type','Spearman');
                end
            end
        end
    else
        %split on second biomarker
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
                    
                    B1 = B_values1(m,Group1{n}); % biomarker for yes group
                    B2 = B_values1(m,Group2{n}); % biomarker for no group
                    [indef,Pvalues(n,m),indef] = ttest2(B1,B2);
                    rho(n,m) = mean(B1)-mean(B2);
                else
                    %split on value
                    Group1{n} = find(B_values2(n,:)>splitValue);% yes
                    Group2{n} = find(B_values2(n,:)<splitValue);% no
                    
                    B1 = B_values1(m,Group1{n}); % biomarker for yes group
                    B2 = B_values1(m,Group2{n}); % biomarker for no group
                    [indef,Pvalues(n,m),indef] = ttest2(B1,B2);
                    rho(n,m) = mean(B1)-mean(B2);
                end
            end
        end
    end
    
    assignin('base','Pvalues',Pvalues)
    assignin('base','rho',rho)
    
    Pvalues(isinf(log10(Pvalues))) = 10^-2.6;
    
    if (display_ind == 1)
        h2 = figure('Visible','on','numbertitle','off','Name','Biomarkers Comparison','position',[10 80 1700 500]);
        %display as grid system
        
        subplot(1,1,1)
        %--- bar plot (map of pvalues)
        minPValue = -2.6;% Plot log10(P-Values) to trick colour bar
        maxPValue = -0;
        hh=uicontextmenu;
        hh2 = uicontextmenu;
        bh=bar3(log10(Pvalues'));
        for k=1:length(bh)
            zdata = get(bh(k),'Zdata');
            set(bh(k),'cdata',zdata);
        end
        colorbar('off')
        coolWarm = load('nbt_CoolWarm.mat','coolWarm');
        coolWarm = coolWarm.coolWarm;
        colormap(coolWarm);
        cbh = colorbar('SouthOutside');
        caxis([minPValue maxPValue])
        set(cbh,'XTick',[-2.6 -1.3010 -1 0])
        set(cbh,'XTicklabel',[0.01 0.05 0.1 1]) %(log scale)
        set(get(cbh,'title'),'String','P-values');
        axis tight
        set(gca,'xticklabel','')
        if strcmp(bioms_name1,'rsq.Answers') == 0
            for j = 1: size(B_values1,1)
                if strcmp(regs_or_chans_name,'Regions')
                    umenu = text(size(B_values2,1)+1,j,[num2str(j) '. ' G(1).chansregs.listregdata(j).reg.name],'horizontalalignment','left','fontsize',8,'interpreter','none','rotation',-30);
                else
                    umenu = text(size(B_values2,1)+1,j,['Channel '  num2str(j)],'horizontalalignment','left','fontsize',8,'rotation',-30);
                end
            end
            
            
            set(gca,'yTick',[],'yticklabel',[],'fontsize',8)
        else
            for j = 1: size(B_values1,1)
                umenu = text(size(B_values2,1)+1,j,[rsq.Questions{j} ' ' num2str(j)],'horizontalalignment','left','fontsize',8,'rotation',-90);
                set(umenu,'uicontextmenu',hh);
            end
            set(gca,'yTick',[],'yticklabel',[],'fontsize',8)
        end
        
        view(-90,-90)
        
        %     pos = get(gca,'Position')
        pos = get(cbh,'Position');
        set(cbh,'Position',[pos(1)-0.5*pos(1) pos(2)-0.6*pos(2) 0.1 0.03])
        for j = 1: size(B_values2,1)
            if strcmp(bioms_name2,'rsq.Answers')
                umenu = text(j,-20,[rsq.Questions{j} ' ' num2str(j)],'horizontalalignment','left','fontsize',8);
                set(umenu,'uicontextmenu',hh);
            else
                if strcmp(regs_or_chans_name,'Regions')
                    umenu = text(j,-5,[num2str(j) '. ' G(1).chansregs.listregdata(j).reg.name],'horizontalalignment','left','fontsize',8,'interpreter','none');
                else
                    umenu = text(j,-12,['Channel '  num2str(j)],'horizontalalignment','left','fontsize',8);
                end
            end
        end
        title(['P-values of the correlation between ',regexprep(bioms2,'_',' '), ' and ', regexprep(bioms1,'_',' ')],'fontweight','bold','fontsize',12)
        title(['Correlation between groups difference of ',regexprep(bioms2,'_',' '), ' and difference ', regexprep(bioms1,'_',' ')],'fontweight','bold','fontsize',12)
        set(bh,'uicontextmenu',hh2);
        
        uimenu(hh,'label','Correlation topoplot','callback',{@nbt_compareBiomarkersPlotTopos,G,B_values1,B_values2,bioms1,bioms2,1,Pvalues',rho',length(group_ind),splitType,splitValue,regs_or_chans_name});
        uimenu(hh2,'label','plottest','callback',{@nbt_compareBiomarkersPlotChansComp,B_values1,B_values2,bioms1,bioms2,1,length(group_ind),splitType,splitValue,Pvalues});
        
        
        
    else
        %display as topoplots
        
        if (strcmp(bioms_name2,'rsq.Answers') == 1)
            %sbp = ceil(size(B_values2,1)^0.5);
            sbp = 4;
            mxlim = max([abs(max(max(rho))) abs(min(min(rho)))]);
            noFigs = ceil(size(B_values2,1)/16);
            for k = 1:noFigs
                figure;
                coolWarm = load('nbt_CoolWarm.mat','coolWarm');
                coolWarm = coolWarm.coolWarm;
                colormap(coolWarm);
                for i = 1+((k-1)*16):min(k*16,size(B_values2,1))
                    idd = i-((k-1)*16);
                    % subplot(sbp,sbp,idd);
                    if strcmp(regs_or_chans_name,'Regions')
                        
                        if nnz(isnan(Pvalues(i,:))) > 0
                            xidd =  mod(idd-1,4);
                            yidd = 4-(1 +floor((idd-1)/4));
                            
                            axes;
                            
                            set(gca,'position',[0.05+ (0.25 * xidd) 0.05 + (0.25 * yidd) 0.2 0.1]);
                            set(gca,'xtick',[]);
                            set(gca,'ytick',[]);
                            text(0.25,0.5,'No members in one of the groups','fontweight','bold');
                            origTitle = get(get(gca,'title'),'position');
                        else
                            
                            axes;
                            set(gca,'clim',[-2.6 0]);
                            nbt_plot_subregions_hack(log10(Pvalues(i,:)),-2.6,0);
                            xidd =  mod(idd-1,4);
                            yidd = 4-(1 +floor((idd-1)/4));
                            
                            cbh = colorbar('WestOutside');
                            oldCBH = get(cbh,'position');
                            oldCBH(1) = 0.05 + (0.25 * xidd) + 0.01;
                            oldCBH(2) = 0.05 + (0.25 * yidd);
                            oldCBH(3) = 0.0086;
                            oldCBH(4) = 0.1009;
                            
                            set(cbh,'position',oldCBH);
                            set(cbh,'yTick',[-2 -1.3010 -1 0])
                            set(cbh,'yTicklabel',[0.01 0.05 0.1 1])
                            set(gca,'position',[0.05 + (0.25 * xidd) 0.05 + (0.25 * yidd) 0.1 0.1]);
                            
                            
                            axes;
                            
                            nbt_plot_subregions_hack(rho(i,:),-mxlim,mxlim);
                            set(gca,'clim',[-mxlim mxlim]);
                            cbh = colorbar('EastOutside');
                            oldCBH = get(cbh,'position');
                            oldCBH(1) = 0.05+ (0.25 * xidd)+0.15;
                            oldCBH(2) = 0.05 + (0.25 * yidd);
                            oldCBH(3) = 0.0086;
                            oldCBH(4) = 0.1009;
                            set(cbh,'position',oldCBH);
                            set(gca,'position',[ 0.05+(0.25 * xidd)+0.07 0.05 + (0.25 * yidd) 0.1 0.1]);
                            origTitle = get(get(gca,'title'),'position');
                            if (test_ind ==1)
                                text(-0.2,-1.2514,origTitle(3),'rho');
                            else
                                text(-0.2,-1.2514,origTitle(3),'\delta mean');
                            end
                            text(-3,-1.2514,origTitle(3),'P-value');
                            origTitle(1) = -2.7;
                            
                            
                            
                        end
                        
                        
                        if (length(rsq.Questions{i}) > 20)
                            bs = rsq.Questions{i};
                            [ab, cd] = strtok(bs(15:end));
                            ad = length(ab);
                            
                            if nnz(isspace(cd)) == length(cd)
                                title([num2str(i) '. ' bs(1:14 + ad)],'FontWeight','Bold');
                            else
                                title({[num2str(i) '. ' bs(1:14 + ad)], cd},'FontWeight','Bold');
                            end
                        else
                            title([num2str(i) '. ' rsq.Questions{i}],'FontWeight','Bold');
                        end
                        
                        
                        
                        set(get(gca,'title'),'position',origTitle);
                        
                    else
                        
                        
                        if nnz(isnan(Pvalues(i,:))) > 0
                            xidd =  mod(idd-1,4);
                            yidd = 4-(1 +floor((idd-1)/4));
                            
                            axes;
                            
                            set(gca,'position',[0.05+ (0.25 * xidd) 0.05 + (0.25 * yidd) 0.2 0.1]);
                            set(gca,'xtick',[]);
                            set(gca,'ytick',[]);
                            text(0.25,0.5,'No members in one of the groups','fontweight','bold');
                            
                        else
                            axes;
                            topoplot(log10(Pvalues(i,:)),G(1).chansregs.chanloc);
                            xidd =  mod(idd-1,4);
                            yidd = 4-(1 +floor((idd-1)/4));
                            set(gca,'clim',[-2.6 0]);
                            
                            cbh = colorbar('WestOutside');
                            oldCBH = get(cbh,'position');
                            oldCBH(1) = 0.05+ (0.25 * xidd) + 0.01;
                            oldCBH(2) = 0.05 + (0.25 * yidd);
                            oldCBH(3) = 0.0086;
                            oldCBH(4) = 0.1009;
                            set(cbh,'position',oldCBH);
                            set(cbh,'yTick',[-2 -1.3010 -1 0])
                            set(cbh,'yTicklabel',[0.01 0.05 0.1 1])
                            
                            
                            set(gca,'position',[0.05+ (0.25 * xidd) 0.05 + (0.25 * yidd) 0.1 0.1]);
                            
                            axes;
                            topoplot(rho(i,:),G(1).chansregs.chanloc);
                            set(gca,'clim',[-mxlim mxlim]);
                            cbh = colorbar('EastOutside');
                            oldCBH = get(cbh,'position');
                            oldCBH(1) = 0.05+(0.25 * xidd)+0.15;
                            oldCBH(2) =  0.05 + (0.25 * yidd);
                            oldCBH(3) = 0.0086;
                            oldCBH(4) = 0.1009;
                            set(cbh,'position',oldCBH);
                            
                            set(gca,'position',[ 0.05+(0.25 * xidd)+0.07 0.05 + (0.25 * yidd) 0.1 0.1]);
                            origTitle = get(get(gca,'title'),'position');
                            if (test_ind ==1)
                                text(-0.1,-origTitle(2),origTitle(3),'rho');
                            else
                                text(-0.1,-origTitle(2),origTitle(3),'\delta mean');
                            end
                            text(-1.4,-origTitle(2),origTitle(3),'P-value');
                            origTitle(1) = -0.7;
                            
                            set(get(gca,'title'),'position',origTitle);
                        end
                        if (length(rsq.Questions{i}) > 20)
                            bs = rsq.Questions{i};
                            [ab, cd] = strtok(bs(15:end));
                            ad = length(ab);
                            
                            if nnz(isspace(cd)) == length(cd)
                                title([num2str(i) '. ' bs(1:14 + ad)],'FontWeight','Bold');
                            else
                                title({[num2str(i) '. ' bs(1:14 + ad)], cd},'FontWeight','Bold');
                            end
                        else
                            title([num2str(i) '. ' rsq.Questions{i}],'FontWeight','Bold');
                        end
                        origTitle = get(get(gca,'title'),'position');
                        
                        origTitle(1) = -0.7;
                        
                        set(get(gca,'title'),'position',origTitle);
                        
                    end
                    
                end
                
                
            end
            
        else
            
            if (strcmp(bioms_name1,'rsq.Answers') == 1)
                
                
                
                sbp = 4;
                mxlim = max([abs(max(max(rho))) abs(min(min(rho)))]);
                noFigs = ceil(size(B_values1,1)/16);
                for k = 1:noFigs
                    figure;
                    coolWarm = load('nbt_CoolWarm.mat','coolWarm');
                    coolWarm = coolWarm.coolWarm;
                    colormap(coolWarm);
                    for i = 1+((k-1)*16):min(k*16,size(B_values1,1))
                        try
                            idd = i-((k-1)*16);
                            % subplot(sbp,sbp,idd);
                            if strcmp(regs_or_chans_name,'Regions')
                                axes;
                                set(gca,'clim',[-2.6 0]);
                                nbt_plot_subregions_hack(log10(Pvalues(:,i)),-2.6,0);
                                xidd =  mod(idd-1,4);
                                yidd = 4-(1 +floor((idd-1)/4));
                                
                                cbh = colorbar('WestOutside');
                                oldCBH = get(cbh,'position');
                                oldCBH(1) = 0.05 + (0.25 * xidd) + 0.01;
                                oldCBH(2) = 0.05 + (0.25 * yidd);
                                oldCBH(3) = 0.0086;
                                oldCBH(4) = 0.1009;
                                
                                set(cbh,'position',oldCBH);
                                set(cbh,'yTick',[-2 -1.3010 -1 0])
                                set(cbh,'yTicklabel',[0.01 0.05 0.1 1])
                                set(gca,'position',[0.05 + (0.25 * xidd) 0.05 + (0.25 * yidd) 0.1 0.1]);
                                
                                
                                axes;
                                
                                nbt_plot_subregions_hack(rho(:,i),-mxlim,mxlim);
                                set(gca,'clim',[-mxlim mxlim]);
                                cbh = colorbar('EastOutside');
                                oldCBH = get(cbh,'position');
                                oldCBH(1) = 0.05+ (0.25 * xidd)+0.15;
                                oldCBH(2) = 0.05 + (0.25 * yidd);
                                oldCBH(3) = 0.0086;
                                oldCBH(4) = 0.1009;
                                set(cbh,'position',oldCBH);
                                set(gca,'position',[ 0.05+(0.25 * xidd)+0.07 0.05 + (0.25 * yidd) 0.1 0.1]);
                                
                                if (length(rsq.Questions{i}) > 20)
                                    bs = rsq.Questions{i};
                                    [ab, cd] = strtok(bs(15:end));
                                    ad = length(ab);
                                    
                                    if nnz(isspace(cd)) == length(cd)
                                        title([num2str(i) '. ' bs(1:14 + ad)],'FontWeight','Bold');
                                    else
                                        title({[num2str(i) '. ' bs(1:14 + ad)], cd},'FontWeight','Bold');
                                    end
                                else
                                    title([num2str(i) '. ' rsq.Questions{i}],'FontWeight','Bold');
                                end
                                
                                origTitle = get(get(gca,'title'),'position');
                                if (test_ind ==1)
                                    text(-0.2,-1.2514,origTitle(3),'rho');
                                else
                                    text(-0.2,-1.2514,origTitle(3),'\delta mean');
                                end
                                text(-3,-1.2514,origTitle(3),'P-value');
                                origTitle(1) = -2.7;
                                
                                set(get(gca,'title'),'position',origTitle);
                                
                            else
                                axes;
                                topoplot(log10(Pvalues(:,i)),G(1).chansregs.chanloc);
                                xidd =  mod(idd-1,4);
                                yidd = 4-(1 +floor((idd-1)/4));
                                set(gca,'clim',[-2.6 0]);
                                
                                cbh = colorbar('WestOutside');
                                oldCBH = get(cbh,'position');
                                oldCBH(1) = 0.05+ (0.25 * xidd) + 0.01;
                                oldCBH(2) = 0.05 + (0.25 * yidd);
                                oldCBH(3) = 0.0086;
                                oldCBH(4) = 0.1009;
                                set(cbh,'position',oldCBH);
                                set(cbh,'yTick',[-2 -1.3010 -1 0])
                                set(cbh,'yTicklabel',[0.01 0.05 0.1 1])
                                
                                
                                set(gca,'position',[0.05+ (0.25 * xidd) 0.05 + (0.25 * yidd) 0.1 0.1]);
                                
                                axes;
                                topoplot(rho(:,i),G(1).chansregs.chanloc);
                                set(gca,'clim',[-mxlim mxlim]);
                                cbh = colorbar('EastOutside');
                                oldCBH = get(cbh,'position');
                                oldCBH(1) = 0.05+(0.25 * xidd)+0.15;
                                oldCBH(2) =  0.05 + (0.25 * yidd);
                                oldCBH(3) = 0.0086;
                                oldCBH(4) = 0.1009;
                                set(cbh,'position',oldCBH);
                                
                                set(gca,'position',[ 0.05+(0.25 * xidd)+0.07 0.05 + (0.25 * yidd) 0.1 0.1]);
                                if (length(rsq.Questions{i}) > 20)
                                    bs = rsq.Questions{i};
                                    [ab, cd] = strtok(bs(15:end));
                                    ad = length(ab);
                                    
                                    if nnz(isspace(cd)) == length(cd)
                                        title([num2str(i) '. ' bs(1:14 + ad)],'FontWeight','Bold');
                                    else
                                        title({[num2str(i) '. ' bs(1:14 + ad)], cd},'FontWeight','Bold');
                                    end
                                else
                                    title([num2str(i) '. ' rsq.Questions{i}],'FontWeight','Bold');
                                end
                                
                                origTitle = get(get(gca,'title'),'position');
                                if (test_ind ==1)
                                    text(-0.1,-origTitle(2),origTitle(3),'rho');
                                else
                                    text(-0.1,-origTitle(2),origTitle(3),'\delta mean');
                                end
                                text(-1.4,-origTitle(2),origTitle(3),'P-value');
                                origTitle(1) = -0.7;
                                
                                set(get(gca,'title'),'position',origTitle);
                            end
                        catch
                        end
                    end
                    
                    
                end
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
            else
                disp('ERROR:  THIS CODE IS NOT WORKING YET');
                %             nbt_plotInsetTopo(([G(1).chansregs.chanloc.X;G(1).chansregs.chanloc.Y;G(1).chansregs.chanloc.Z]'),log10(Pvalues)',[-2.6 0]);
                %             coolWarm = load('nbt_CoolWarm.mat','coolWarm');
                %             coolWarm = coolWarm.coolWarm;
                %             colormap(coolWarm);
                %
                %
                %             mxlim = max([abs(max(max(rho))) abs(min(min(rho)))]);
                %
                %             nbt_plotInsetTopo(([G(1).chansregs.chanloc.X;G(1).chansregs.chanloc.Y;G(1).chansregs.chanloc.Z]'),rho',[-mxlim mxlim]);
                %             coolWarm = load('nbt_CoolWarm.mat','coolWarm');
                %             coolWarm = coolWarm.coolWarm;
                %             colormap(coolWarm);
                
                
            end
            
            
            
            
        end
        
        
        
    end
    
    
    
    
    
    %%See what test if being run
    
    
    
    
    
    
end
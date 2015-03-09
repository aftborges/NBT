function nbt_splitPlotGroups(d1,d2,ListSplit,ListSplitValue,ListBiom1, ListBiom2, ListRegion, ListGroup,G)
global Factors
global Questionnaire
splt = get(ListSplit,'Value');
spltVal = str2num(get(ListSplitValue,'String'));


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



%%  If more than one group then has to be paired and make sure that
%  both have the same subjects.
if length(group_ind) == 1
    Group = G(group_ind);
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
    
elseif length(group_ind) == 2
    Group1 = G(group_ind(1));
    Group2 = G(group_ind(2));
    n_files1 = length(Group1.fileslist);
    n_files2 = length(Group2.fileslist);
    bioms1 = bioms_name1;
    bioms2 = bioms_name2;
    
    nameG2 = Group2.fileslist.group_name;
    path2 = Group2.fileslist.path;
    n_files2 = length(Group2.fileslist);
    nameG1 = Group1.fileslist.group_name;
    path1 = Group1.fileslist.path;
    n_files1 = length(Group1.fileslist);
    % load biomarkers
    for j = 1:n_files1 % subject
        for l = 1:length(bioms2) % biomarker
            namefile = Group1.fileslist(j).name;
            biomarker1 = bioms1{1};
            biomarker2 = bioms2{l};
            [B1_values1(:,j),Sub,Proj,unit1{j}]=nbt_load_analysis(path1,namefile,biomarker1,@nbt_get_biomarker,[],[],[]);
            [B1_values2(:,j,l),Sub,Proj,unit2{j,l}]=nbt_load_analysis(path1,namefile,biomarker2,@nbt_get_biomarker,[],[],[]);
            
        end
    end
    for j = 1:n_files2 % subject
        for l = 1:length(bioms2) % biomarker
            namefile = Group2.fileslist(j).name;
            biomarker1 = bioms1{1};
            biomarker2 = bioms2{l};
            [B2_values1(:,j),Sub,Proj,unit1{j}]=nbt_load_analysis(path2,namefile,biomarker1,@nbt_get_biomarker,[],[],[]);
            [B2_values2(:,j,l),Sub,Proj,unit2{j,l}]=nbt_load_analysis(path2,namefile,biomarker2,@nbt_get_biomarker,[],[],[]);
        end
    end
    
    B_values1 = B1_values1-B2_values1;% questionnaire
    B_values2 = B1_values2-B2_values2;%i.e. amplitude
    
end

if strcmp(regs_or_chans_name,'Regions')
    if isempty(strfind(bioms_name1,'.Answers'))
        regions = G(group_ind(1)).chansregs.listregdata;
        for j = 1:size(B_values1,2) % subject
            
            B1 = B_values1(:,j);
            B_gebruik1(:,j) = nbt_compare_getRegions(B1,regions);
        end
        clear B_values1;
        B_values1 = B_gebruik1;
    end
    if isempty(strfind(bioms_name2,'.Answers'))
        regions = G(group_ind(1)).chansregs.listregdata;
        for j = 1:size(B_values2,2) % subject
            B2 = B_values2(:,j);
            B_gebruik2(:,j) = nbt_compare_getRegions(B2,regions);
        end
        clear B_values2;
        B_values2 = B_gebruik2;
    end
end

vals = B_values2;
nm = bioms_name2;

%sbp = ceil(size(vals,1)^0.5);
sbp = 4; %now showing 16 graphs per figure;
noSubj = size(vals,2);
noFigs = ceil(size(vals,1)/16);

for k = 1:noFigs
    figure
    for i = 1+(k-1)*16:min(k*16,size(vals,1));
        %plot cdf
        %indicate high and low groups
      hSplit = subplot(sbp,sbp,i-(k-1)*16);
      
        plot(sort(vals(i,:)),1/noSubj:1/noSubj:1);
        
        if ~isempty(strfind(nm,'.Answers'))
            if cellfun(@isempty,(strfind(nm,'Factors')))
                varQuest = Questionnaire;
            else
                varQuest = Factors;
            end
            if (length(varQuest.Questions{i}) > 20)
                bs = varQuest.Questions{i};
                [ab, cd] = strtok(bs(15:end));
                ad = length(ab);
                
                if nnz(isspace(cd)) == length(cd)
                    title([num2str(i) '. ' bs(1:14 + ad)],'FontWeight','Bold');
                else
                    title({[num2str(i) '. ' bs(1:14 + ad)], cd},'FontWeight','Bold');
                end
            else
                title([num2str(i) '. ' varQuest.Questions{i}],'FontWeight','Bold');
            end
            if length(group_ind) ==1
                xlabel([strtok(nm,'.') ' score']);
            else
                xlabel(['Change in ' strtok(nm,'.') ' score']);
            end
            ylabel('CDF');
            
        else
            if strcmp(regs_or_chans_name,'Regions')
                title([num2str(i) '. ' G(1).chansregs.listregdata(i).reg.name],'interpreter','none');
            else
            title(['Channel ' num2str(i)]);
            end
            if length(group_ind) ==1
                xlabel([strtok(bioms_name2,'.') ' Value'],'interpreter','none');
            else
                xlabel(['Change in ' strtok(bioms_name2,'.') ' Value'],'interpreter','none');
            end
            ylabel('CDF');
        end
        axis([min(min(vals)) max(max(vals)) 0 1]);     
        if splt == 1
            
           hlow  = refline(0,spltVal/100);
           hhigh =  refline(0,1-(spltVal/100));
           set(hhigh,'Color','r');
           set(hlow,'Color','g');
            
        else
            vs = sort(vals(i,:));
            if (nnz(vs>spltVal) == 0 || nnz(vs<spltVal) == 0)
                pos = get(gca,'Position');
                midY = 0.5;
                midX = 0.5;
                text(midX,midY,['All values either above ';'or below the split value']);
                axis off;
            else
            mx = min(find(vs>spltVal));
            mn = max(find(vs<spltVal));
            hhigh = refline(0,mx/noSubj);
            hlow = refline(0,mn/noSubj);
            set(hhigh,'Color','r');
            set(hlow,'Color','g');
                
            end
            
            
        end
        set(hSplit,'ButtonDownFcn',@nbt_plotLargeSplitGraph);
    end
    
   
   
end

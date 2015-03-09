function [B_values_cell,Sub,Proj,unit] = nbt_checkif_groupdiff(Group,G,n_files,bioms_name,path);

if  ~isempty(Group.group_difference)% check if the group comes froma group difference 
   diffname = Group.group_difference;
   findminus = findstr(diffname,'-');
   g_diff1 = strtrim(diffname(1:findminus-1)); 
   g_diff2 = strtrim(diffname(findminus+1:end)); 
   for i = 1:length(G)
       if strcmp(strtrim(G(i).fileslist(1).group_name),g_diff1)
           firstG = i;
           break
       end
   end
   for i = 1:length(G)
       if strcmp(strtrim(G(i).fileslist(1).group_name),g_diff2)
           secondG = i;
           break
       end
   end
   firstG = G(firstG);
   secondG = G(secondG);
if n_files ~= length(firstG.fileslist)
    n_files =length(firstG.fileslist);
end
   for j = 1:n_files % subjects
    disp(j)
    for l = 1:length(bioms_name) % biomarker
        namefilefirstG = firstG.fileslist(j).name;
        [B_values_cellfirstG{j,l},Sub,Proj,unit{j,l}]=nbt_load_analysis(path,namefilefirstG,bioms_name{l},@nbt_get_biomarker,[],[],[]);
        namefilesecondG = secondG.fileslist(j).name;
        [B_values_cellsecondG{j,l},Sub,Proj,unit{j,l}]=nbt_load_analysis(path,namefilesecondG,bioms_name{l},@nbt_get_biomarker,[],[],[]);
        B_values_cell{j,l} = B_values_cellfirstG{j,l}-B_values_cellsecondG{j,l};
    end
  end
           
           
else
    for j = 1:n_files % subjects
        disp(j)
        for l = 1:length(bioms_name) % biomarker
            namefile = Group.fileslist(j).name;
            [B_values_cell{j,l},Sub,Proj,unit{j,l}]=nbt_load_analysis(path,namefile,bioms_name{l},@nbt_get_biomarker,[],[],[]);
        end
    end
end

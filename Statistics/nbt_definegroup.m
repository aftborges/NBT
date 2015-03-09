% nbt_definegroup - Open Interface for defining a statistical group,
%                   it generates a new directory struct called SelectedFiles
%                   which contains only the files selected for the statistics
%
% Usage:
%   nbt_definegroup(d)
%
% Inputs:
%   d is a directory struct obtained as follow d = dir(path); with path
%   indicating the location of NBT files
%
% Outputs:
%
% Example:
%
% References:
%
% See also:
%

%------------------------------------------------------------------------------------
% Originally created by Giuseppina Schiavone (2012), see NBT website (http://www.nbtwiki.net) for current email address
%------------------------------------------------------------------------------------
%
% ChangeLog - see version control log at NBT website for details.
%
% Copyright (C) <year>  <Main Author>  (Neuronal Oscillations and Cognition group,
% Department of Integrative Neurophysiology, Center for Neurogenomics and Cognitive Research,
% Neuroscience Campus Amsterdam, VU University Amsterdam)
%
% Part of the Neurophysiological Biomarker Toolbox (NBT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% See Readme.txt for additional copyright information.
% -------------------------------------------------------------------------


function nbt_definegroup(varargin)
P = varargin;
nargs = length(P);


con = 1;
sub = 1;
pro = 1;
date = 1;
gender = 1;
age = 1;
readconditions = '';
readproject = '';
readsubject = '';
readdate = '';
readgender = {''};
readage = {''};

if(nbt_determineNBTelementState) %use NBTelement
    %load NBTelementBase from base
    readconditions = evalin('base', 'Condition.Data');
    SubjList = evalin('base', 'Subject.Data');
    for mm=1:length(SubjList);
        readsubject{mm} = num2str(SubjList(mm));
    end
else %then we read the information from the analysis files
    if nargs<1
        [path]=uigetdir([],'Select folder with NBT Signals');
        d = dir(path);
    else
        d = P{1};
    end
    
    %--- scan files in the folder
    %--- for files copied from a mac
    startindex = 0;
    for i = 1:length(d)
        if  d(i).isdir || strcmp(d(i).name(1),'.') || strcmp(d(i).name(1:2),'..') || strcmp(d(i).name(1:2),'._')
            startindex = i+1;
        end
    end
    %---
    
    for i = startindex:length(d)
        if isempty(findstr(d(i).name,'analysis')) && ~isempty(findstr(d(i).name,'info')) && ~isempty(findstr(d(i).name(end-3:end),'.mat')) && isempty(findstr(d(i).name,'statistics'))
            index = findstr(d(i).name,'.');
            index2 = findstr(d(i).name,'_');
            % read gender and age
            % load info file
            Loaded = load([path '/' d(i).name]);
            
            Infofields = fieldnames(Loaded);
            Firstfield = Infofields{1};
            clear Loaded Infofields;
            
            SignalInfo = load([path '/' d(i).name],Firstfield);
            
            SignalInfo = eval(strcat('SignalInfo.',Firstfield));
            disp('Please wait: NBT is checking the files into your folder...')
            if i == startindex
                readconditions{con} =  d(i).name(index(3)+1:index2-1);
                readproject{pro} =  d(i).name(1:index(1)-1);
                readsubject{sub} = d(i).name(index(1)+1:index(2)-1);
                readdate{date} = d(i).name(index(2)+1:index(3)-1);
                if ~isempty(SignalInfo.subject_gender)
                    readgender{gender} = SignalInfo.subject_gender;
                else
                    readgender{gender} = [];
                end
                if ~isempty(SignalInfo.subject_gender)
                    readage{age} = SignalInfo.subject_age;
                else
                    readage{age} = [];
                end
                
                
                con = con+1;
                pro = pro+1;
                sub = sub+1;
                date = date+1;
                gender = gender +1;
                age = age+1;
            else
                if ~strcmp(readconditions,d(i).name(index(3)+1:index2-1))
                    readconditions{con} = d(i).name(index(3)+1:index2-1);
                    con = con+1;
                end
                if ~strcmp(readproject,d(i).name(1:index(1)-1))
                    readproject{pro} = d(i).name(1:index(1)-1);
                    pro = pro+1;
                end
                if ~strcmp(readsubject,d(i).name(index(1)+1:index(2)-1))
                    readsubject{sub} = d(i).name(index(1)+1:index(2)-1);
                    sub = sub+1;
                end
                if ~strcmp(readdate,d(i).name(index(2)+1:index(3)-1))
                    readdate{date} = d(i).name(index(2)+1:index(3)-1);
                    date = date+1;
                end
                if ~(length(readgender) == 2)
                    if ~strcmp(num2str(readgender{end}),num2str(SignalInfo.subject_gender))
                        
                        if ~isempty(SignalInfo.subject_gender)
                            readgender{gender} = SignalInfo.subject_gender;
                        else
                            readgender{gender} = [];
                        end
                        gender = gender+1;
                    end
                end
                if ~strcmp(num2str(readage{end}),num2str(SignalInfo.subject_age))
                    if ~isempty(SignalInfo.subject_age)
                        readage{age} = SignalInfo.subject_age;
                    else
                        readage{age} =[];
                    end
                    age = age+1;
                end
            end
        end
    end
    clear con sub date pro age gender
end
if nargs==1;
    selection.con = P{2};
    selection.sub = readsubject;
    selection.date = readdate;
    selection.pro = readproject;
    selection.age = readage;
    selection.gender = readgender;
elseif nargs ==2;
    selection.con = P{2};
    selection.sub = P{3};
    selection.date =  readdate;
    selection.pro = readproject;
    selection.age = readage;
    selection.gender = readgender;
elseif nargs == 3;
    selection.con = P{2};
    selection.sub = P{3};
    selection.date = P{4};
    selection.pro = readproject;
    selection.age = readage;
    selection.gender = readgender;
elseif nargs == 4;
    selection.con = P{2};
    selection.sub = P{3};
    selection.date = P{5};
    selection.pro = P{6};
    selection.age = readage;
    selection.gender = readgender;
end
%--- interface
if nargs == 0
    scrsz = get(0,'ScreenSize');
    GroupSelection = figure('Units','pixels', 'name','Define Group' ,'numbertitle','off','Position',[390.0000  456.7500  1000  300], ...
        'MenuBar','none','NextPlot','new','Resize','off');
    GroupSelection=nbt_movegui(GroupSelection);

    g = gcf;
    Col = get(g,'Color');
    
    listBox1 = uicontrol(GroupSelection,'Style','listbox','Units','characters',...
        'Position',[4 1 15 17],...
        'BackgroundColor','white',...
        'Max',10,'Min',1, 'String', readconditions,'Value',[]);
    get(listBox1,'Position');
    text_ui1= uicontrol(GroupSelection,'Style','text','Position',[20 270 700 20],...
        'string','Condition            Subject             Date                 Project               Gender             Age   ','fontsize',12, 'HorizontalAlignment','Left');
    listBox2 = uicontrol(GroupSelection,'Style','listbox','Units','characters',...
        'Position',[24 1 15 17],...
        'BackgroundColor','white',...
        'Max',10,'Min',1, 'String', readsubject,'Value',[]);
%     text_ui2= uicontrol(GroupSelection,'Style','text','Position',[160 270 100 20],'string',' Subject ','fontsize',12);
    listBox3 = uicontrol(GroupSelection,'Style','listbox','Units','characters',...
        'Position',[44 1 15 17],...
        'BackgroundColor','white',...
        'Max',10,'Min',1, 'String', readdate,'Value',[]);
%     text_ui3= uicontrol(GroupSelection,'Style','text','Position',[260 270 100 20],'string','  Date   ','fontsize',12);
    listBox4 = uicontrol(GroupSelection,'Style','listbox','Units','characters',...
        'Position',[64 1 15 17],...
        'BackgroundColor','white',...
        'Max',10,'Min',1, 'String', readproject,'Value',[]);
%     text_ui4= uicontrol(GroupSelection,'Style','text','Position',[360 270 100 20],'string',' Project ','fontsize',12);
    listBox5 = uicontrol(GroupSelection,'Style','listbox','Units','characters',...
        'Position',[84 1 15 17],...
        'BackgroundColor','white',...
        'Max',10,'Min',1, 'String', readgender,'Value',[]);
%     text_ui5= uicontrol(GroupSelection,'Style','text','Position',[460 270 100 20],'string',' Gender  ','fontsize',12);
    % sort age
    if ~isempty(cell2mat(readage))
    countage = 1;
    if(~strcmp(readage{1,1},''))
    for anni = 1:length(readage) 
        if ischar(readage{anni})
        eta(countage) = str2num(readage{anni});  
        else
            if isempty((readage{anni}))
                eta(countage) = 99;  
            else
                eta(countage) = (readage{anni});  
            end
        end
        countage = countage + 1;     
    end
    eta = sort(eta);
    eta = unique(eta);
    clear readage
    for anni = 1:length(eta)
        readage{anni} = eta(anni);
    end
    end
    end
    
    %
    listBox6 = uicontrol(GroupSelection,'Style','listbox','Units','characters',...
        'Position',[104 1 15 17],...
        'BackgroundColor','white',...
        'Max',10,'Min',1, 'String', readage,'Value',[]);
%     text_ui6= uicontrol(GroupSelection,'Style','text','Position',[560 270 100 20],'string','   Age   ','fontsize',12);
    
<<<<<<< HEAD
    text_ui7= uicontrol(GroupSelection,'Style','text','Position',[700 70 200 20],'string','Write a name for the Group','fontsize',10);
    text_ui8= uicontrol(GroupSelection,'Style','edit','Position',[700 50 200 20],'string','','fontsize',10);
    
    plotButton = uicontrol(GroupSelection,'Style','pushbutton','Units','pixels','Position',[920 30 40 40], 'String','OK','callback', @groupdefinition);
=======
    text_ui7= uicontrol(GroupSelection,'Style','text','Position',[760 270 200 20],'string','Insert a name for the Group','fontsize',10);
    text_ui8= uicontrol(GroupSelection,'Style','edit','Position',[760 230 200 20],'string','','fontsize',10);
    
    plotButton = uicontrol(GroupSelection,'Style','pushbutton','Units','characters','Position',[140 14 20 2], 'String','OK','callback', @groupdefinition);
    
>>>>>>> FETCH_HEAD
else
    
    startindex = 0;
    for i = 1:length(d)
        if  d(i).isdir || strcmp(d(i).name(1),'.') || strcmp(d(i).name(1:2),'..') || strcmp(d(i).name(1:2),'._')
            startindex = i+1;
        end
    end
    k = 1;
    for i = startindex:length(d)
        if ~isempty(findstr('analysis',d(i).name))
            analysis_files(k) = d(i);
            k = k +1;
        end
    end
    g =1;
    for i = 1:length(selection.con)
        for j = 1:length(selection.sub)
            for k = 1:length(selection.date)
                for h = 1:length(selection.pro)
                    SelFiles{g} = strcat(selection.pro(h),'.',selection.sub(j),'.',selection.date(k),'.',selection.con(i),'_analysis.mat');
                    g = g+1;
                end
            end
        end
    end
    
    k = 1;
    for i = startindex:length(d)
        % read gender and age
        for j = 1:length(SelFiles)
            if strcmp(d(i).name,cell2mat(SelFiles{j}))
                Loaded = load([path '/' d(i).name(1:end-13) '_info.mat']);
                
                Infofields = fieldnames(Loaded);
                Firstfield = Infofields{1};
                clear Loaded Infofields;
                
                SignalInfo = load([path '/' d(i).name(1:end-13) '_info.mat'],Firstfield);
                
                SignalInfo = eval(strcat('SignalInfo.',Firstfield));
                %-check gender
                for m = 1:length(selection.gender);
                    if ~ischar(SignalInfo.subject_gender)
                        SignalInfo.subject_gender = num2str(SignalInfo.subject_gender);
                    end
                    if ~ischar(SignalInfo.subject_age)
                        SignalInfo.subject_age = num2str(SignalInfo.subject_age);
                    end
                    if isequal(SignalInfo.subject_gender,num2str(selection.gender{m}))
                        for n = 1:length(selection.age)
                            if isequal((SignalInfo.subject_age),num2str(selection.age{n}))
                                %-check age
                                disp('Please wait: NBT is sorting the files...')
                                d(i).path = path;
                                group_name = get(text_ui8,'String');
                                d(i).group_name = group_name;
                                SelectedFiles(k) = d(i); % contains exactly the files to be used for the statistics
                                k = k +1;
                            end
                        end
                    end
                end
            end
        end
    end
    assignin('base','SelectedFiles',SelectedFiles)
    close all
end

% --- callback function - nested function
    function groupdefinition(src,evt)
        set(plotButton,'String', 'Busy...');
        vars1 = get(listBox1,'String');
        var_index1 = get(listBox1,'Value');
        if isempty(var_index1)
            selection.con = vars1;
        else
            if length(vars1) == length(var_index1)
                selection.con = vars1;
            else
                selection.con = vars1(var_index1);
            end
        end
        
        vars2 = get(listBox2,'String');
        var_index2 = get(listBox2,'Value');
        if isempty(var_index2)
            selection.sub = vars2;
        else
            if length(vars2) == length(var_index2)
                selection.sub = vars2;
            else
                selection.sub = vars2(var_index2);
            end
        end
        
        vars3 = get(listBox3,'String');
        var_index3 = get(listBox3,'Value');
        if isempty(var_index3)
            selection.date = vars3;
        else
            if length(vars3) == length(var_index3)
                selection.date = vars3;
            else
                selection.date = vars3(var_index3);
            end
        end
        
        
        vars4 = get(listBox4,'String');
        var_index4 = get(listBox4,'Value');
        if isempty(var_index4)
            selection.pro = vars4;
        else
            if length(vars4) == length(var_index4)
                selection.pro = vars4;
            else
                selection.pro = vars4(var_index4);
            end
        end
        
        vars5 = get(listBox5,'String');
        var_index5 = get(listBox5,'Value');
        if isempty(var_index5)
            selection.gender = vars5;
        else
            if length(vars5) == length(var_index5)
                selection.gender = vars5;
            else
                selection.gender = vars5(var_index5);
            end
        end
        
        vars6 = get(listBox6,'String');
        var_index6 = get(listBox6,'Value');
        if isempty(var_index6)
            selection.age = vars6;
        else
            if length(vars6) == length(var_index6)
                selection.age = vars6;
            else
                selection.age = vars6(var_index6);
            end
        end
        group_name = get(text_ui8,'String');
        selection.group_name = group_name;
        %         selection.con = con;
        %         selection.sub = sub;
        %         selection.date = date;
        %         selection.pro = pro;
        
        
        if(~nbt_determineNBTelementState) %loading information from analysis files
            %--- generate Selected file dir struct
            startindex = 0;
            for i = 1:length(d)
                if  d(i).isdir || strcmp(d(i).name(1),'.') || strcmp(d(i).name(1:2),'..') || strcmp(d(i).name(1:2),'._')
                    startindex = i+1;
                end
            end
            k = 1;
            for i = startindex:length(d)
                if ~isempty(findstr('analysis',d(i).name))
                    analysis_files(k) = d(i);
                    k = k +1;
                end
            end
            g =1;
            for i = 1:length(selection.con)
                for j = 1:length(selection.sub)
                    for k = 1:length(selection.date)
                        for h = 1:length(selection.pro)
                            SelFiles{g} = strcat(selection.pro(h),'.',selection.sub(j),'.',selection.date(k),'.',selection.con(i),'_analysis.mat');
                            g = g+1;
                        end
                    end
                end
            end
            
            
            k = 1;
            for i = startindex:length(d)
                % read gender and age
                for j = 1:length(SelFiles)
                    if strcmp(d(i).name,cell2mat(SelFiles{j}))
                        
                        Loaded = load([path '/' d(i).name(1:end-13) '_info.mat']);
                        
                        Infofields = fieldnames(Loaded);
                        Firstfield = Infofields{1};
                        clear Loaded Infofields;
                        
                        SignalInfo = load([path '/' d(i).name(1:end-13) '_info.mat'],Firstfield);
                        
                        SignalInfo = eval(strcat('SignalInfo.',Firstfield));
                        %-check gender
                        for m = 1:length(selection.gender);
                            if ~ischar(SignalInfo.subject_gender)
                                SignalInfo.subject_gender = num2str(SignalInfo.subject_gender);
                            end
                            if ~ischar(SignalInfo.subject_age)
                                SignalInfo.subject_age = num2str(SignalInfo.subject_age);
                            end
                            if isequal(SignalInfo.subject_gender,(selection.gender{m}))
                                for n = 1:length(selection.age)
                                    if isequal((SignalInfo.subject_age),num2str(selection.age{n}))
                                        %-check age
                                        disp('Please wait: NBT is sorting the files...')
                                        d(i).path = path;
                                        group_name = get(text_ui8,'String');
                                        d(i).group_name = group_name;
                                        SelectedFiles(k) = d(i); % contains exactly the files to be used for the statistics
                                        k = k +1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            assignin('base','SelectedFiles',SelectedFiles)
        else
            assignin('base', 'selection',selection);
        end
        h = get(0,'CurrentFigure');
        close(h)
    end
end

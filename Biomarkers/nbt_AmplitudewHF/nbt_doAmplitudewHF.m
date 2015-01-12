% nbt_doAmplitudeHF - Computes amplitudes classical frequencies and high
% gamma
%

% Usage:
% [amplitude_1_4_Hz amplitude_4_8_Hz amplitude_8_13_Hz amplitude_13_30_Hz amplitude_30_45_Hz ... 
%  amplitude_1_4_Hz_Normalized amplitude_4_8_Hz_Normalized amplitude_8_13_Hz_Normalized  ...
%  amplitude_13_30_Hz_Normalized amplitude_30_45_Hz_Normalized] = nbt_doAmplitude(Signal,SignalInfo);
%

% NOTE: omitting the normalized amplitude version for now
% Inputs:
%
% Signal = NBT Signal matrix
% SignalInfo = NBT Info object
%
% Outputs: several amplitude biomarker objects
%
% This function creates amplitude biomarker objects, where it stores integrated amplitudes in several frequency bands per channel and per sub region.
% Optionally, it plots integrated amplitudes in several frequency bands in
% three different ways: per channel color coded, spatially interpolated,
% and the mean in 6 subregions color coded.

%------------------------------------------------------------------------------------
% Originally created by Rick Jansen (2012), see NBT website (http://www.nbtwiki.net) for current email address
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
% function [amplitude_1_4_Hz amplitude_4_8_Hz amplitude_8_13_Hz amplitude_13_30_Hz amplitude_30_45_Hz ...
%     amplitude_1_4_Hz_Normalized amplitude_4_8_Hz_Normalized amplitude_8_13_Hz_Normalized  ...
%     amplitude_13_30_Hz_Normalized amplitude_30_45_Hz_Normalized,amplitude_1_4_Hz_cmbgrad,amplitude_4_8_Hz_cmbgrad, ...
%     amplitude_8_13_Hz_cmbgrad,amplitude_13_30_Hz_cmbgrad,amplitude_30_45_Hz_cmbgrad,amplitude_1_4_Hz_cmbgradN, ...
%     amplitude_4_8_Hz_cmbgradN,amplitude_8_13_Hz_cmbgradN,amplitude_13_30_Hz_cmbgradN,amplitude_30_45_Hz_cmbgradN] = nbt_doAmplitude(Signal,Info)
function [amplitude_1_4_Hz,amplitude_4_8_Hz,amplitude_8_13_Hz,amplitude_13_30_Hz,...
    amplitude_30_45_Hz,amplitude_55_125_Hz,amplitude_60_90_Hz,amplitude_1_4_Hz_Normalized,amplitude_4_8_Hz_Normalized,amplitude_8_13_Hz_Normalized,amplitude_13_30_Hz_Normalized,amplitude_30_45_Hz_Normalized,amplitude_55_125_Hz_Normalized,amplitude_60_90_Hz_Normalized] = nbt_doAmplitudewHF(Signal,Info)
%--- input checks
% error(nargchk(4,4,nargin))

%%   give information to the user
disp(['Computing Amplitudes for ',Info.file_name])

%%    assigning fields
coolWarm = load('nbt_CoolWarm','coolWarm');
coolWarm = coolWarm.coolWarm;
color_scale = '0-max'; %'min-max'; % '0-max' # set color scale: from zero to max or from min to max
save_pdf_plot=0;
plotting=0;
nfft=2^10; %number of fast fourier transforms, higher this number and the frequency resolution of the spectrum goes up,
number_of_channels=size(Signal,2);

%% remove artifact intervals
Signal = nbt_RemoveIntervals(Signal,Info);

%% determine intervals at which power is calculated

FS=Info.converted_sample_frequency;


interval_Hz(1,:)=[1 4];
interval_Hz(2,:)=[4 8];
interval_Hz(3,:)=[8 13];
interval_Hz(4,:)=[13 30];
interval_Hz(5,:)=[30 45];
interval_Hz(6,:)=[55 125];
interval_Hz(7,:)=[60 90];

[p,f]=pwelch(randn(1,nfft),hamming(nfft),0,nfft,FS);

for i=1:7
    interval{i}=find(f>interval_Hz(i,1)&f<interval_Hz(i,2));
end

nr_interv=length(interval);

%% calculated integrated and normalized power for all channels

integrated=zeros(nr_interv,number_of_channels);
normalized=zeros(nr_interv,number_of_channels);




for i=1:number_of_channels
    
    
    
    [p,f]=pwelch(Signal(:,i),hamming(nfft),0,nfft,FS);
    p=sqrt(p); % transform power to amplitude
    
  
   


     for j=1:nr_interv
        integrated(j,i) = mean(p(interval{j})); %%so provisoriamente p calc power em metade S1 e S4
        
        normalized(j,i) = mean(p(interval{j}))/sum([mean(p([interval{1}]));mean(p([interval{2}]));mean(p([interval{3}]));mean(p([interval{4}]));mean(p([interval{5}]))]);
        
    end
end

%% Set Bad Channels to NaNs
integrated(:,find(Info.BadChannels)) = NaN;
normalized(:,find(Info.BadChannels)) = NaN;

%% Create amplitude objects and assign values
 
 amplitude_1_4_Hz=nbt_amplitude(number_of_channels,integrated(1,:),[],0,'\muV');
 amplitude_4_8_Hz=nbt_amplitude(number_of_channels,integrated(2,:),[],0,'\muV');
 amplitude_8_13_Hz=nbt_amplitude(number_of_channels,integrated(3,:),[],0,'\muV');
 amplitude_13_30_Hz=nbt_amplitude(number_of_channels,integrated(4,:),[],0,'\muV');
 amplitude_30_45_Hz=nbt_amplitude(number_of_channels,integrated(5,:),[],0,'\muV');
 amplitude_55_125_Hz=nbt_amplitude(number_of_channels,integrated(6,:),[],0,'\muV');
 amplitude_60_90_Hz=nbt_amplitude(number_of_channels,integrated(7,:),[],0,'\muV');

amplitude_1_4_Hz.FrequencyRange = interval_Hz(1,:);
amplitude_4_8_Hz.FrequencyRange = interval_Hz(2,:);
amplitude_8_13_Hz.FrequencyRange = interval_Hz(3,:);
amplitude_13_30_Hz.FrequencyRange = interval_Hz(4,:);
amplitude_30_45_Hz.FrequencyRange = interval_Hz(5,:);
amplitude_55_125_Hz.FrequencyRange = interval_Hz(6,:);
amplitude_60_90_Hz.FrequencyRange = interval_Hz(7,:);

amplitude_1_4_Hz_Normalized=nbt_amplitude(number_of_channels,normalized(1,:)*100,[],1,'%');
amplitude_4_8_Hz_Normalized=nbt_amplitude(number_of_channels,normalized(2,:)*100,[],1,'%');
amplitude_8_13_Hz_Normalized=nbt_amplitude(number_of_channels,normalized(3,:)*100,[],1,'%');
amplitude_13_30_Hz_Normalized=nbt_amplitude(number_of_channels,normalized(4,:)*100,[],1,'%');
amplitude_30_45_Hz_Normalized=nbt_amplitude(number_of_channels,normalized(5,:)*100,[],1,'%');
amplitude_55_125_Hz_Normalized=nbt_amplitude(number_of_channels,normalized(6,:)*100,[],1,'%');
amplitude_60_90_Hz_Normalized=nbt_amplitude(number_of_channels,normalized(7,:)*100,[],1,'%');

amplitude_1_4_Hz_Normalized.FrequencyRange = interval_Hz(1,:);
amplitude_4_8_Hz_Normalized.FrequencyRange = interval_Hz(2,:);
amplitude_8_13_Hz_Normalized.FrequencyRange = interval_Hz(3,:);
amplitude_13_30_Hz_Normalized.FrequencyRange = interval_Hz(4,:);
amplitude_30_45_Hz_Normalized.FrequencyRange = interval_Hz(5,:);
amplitude_55_125_Hz_Normalized.FrequencyRange = interval_Hz(6,:);
amplitude_60_90_Hz_Normalized.FrequencyRange = interval_Hz(7,:);




%% in case of 129-channel EEG data, compute values at 6 subregions and check plotting option

if isfield(Info.Interface,'EEG') && number_of_channels==129
    
    % integrated(:,129)=nanmedian(integrated'); % set reference channel to mean of all channels
    % normalized(:,129)=nanmedian(normalized');
    
    for ii=1:nr_interv
        integrated_regions(ii,:)= nbt_plot_subregions(integrated(ii,:),0);
        normalized_regions(ii,:)=nbt_plot_subregions(normalized(ii,:),0);
    end
    normalized = normalized * 100;
    normalized_regions = normalized_regions*100;
    
    
    amplitude_1_4_Hz=nbt_amplitude(number_of_channels,integrated(2,:),integrated_regions(1,:),0,'\muV');
    amplitude_4_8_Hz=nbt_amplitude(number_of_channels,integrated(3,:),integrated_regions(2,:),0,'\muV');
    amplitude_8_13_Hz=nbt_amplitude(number_of_channels,integrated(4,:),integrated_regions(3,:),0,'\muV');
    amplitude_13_30_Hz=nbt_amplitude(number_of_channels,integrated(5,:),integrated_regions(4,:),0,'\muV');
    amplitude_30_45_Hz=nbt_amplitude(number_of_channels,integrated(6,:),integrated_regions(5,:),0,'\muV');
    amplitude_55_125_Hz=nbt_amplitude(number_of_channels,integrated(1,:),integrated_regions(6,:),0,'\muV');
    amplitude_60_90_Hz=nbt_amplitude(number_of_channels,integrated(1,:),integrated_regions(7,:),0,'\muV');
    
    amplitude_1_4_Hz.FrequencyRange = interval_Hz(1,:);
    amplitude_4_8_Hz.FrequencyRange = interval_Hz(2,:);
    amplitude_8_13_Hz.FrequencyRange = interval_Hz(3,:);
    amplitude_13_30_Hz.FrequencyRange = interval_Hz(4,:);
    amplitude_30_45_Hz.FrequencyRange = interval_Hz(5,:);
    amplitude_55_125_Hz.FrequencyRange = interval_Hz(6,:);
    amplitude_60_90_Hz.FrequencyRange = interval_Hz(7,:);
    
    
    
    
    %if plotting
        
        %----------- Plot integrated amplitudes
        
        %--- get figure handle
%         
%         done=0;
%         figHandles = findobj('Type','figure');
%         for i=1:length(figHandles)
%             if strcmp(get(figHandles(i),'Tag'),'amplitudes')
%                 figure(figHandles(i))
%                 done=1;
%             end
%         end
%         if ~done
%             figure()
%             set(gcf,'Tag','amplitudes');
%         end
%         set(gcf,'numbertitle','off');
%         set(gcf,'name','NBT: Integrated amplitudes');
%         clf
%         
%         %----plot topoplots per frequency band
%         
%         colormap(coolWarm);
%         for i=1:nr_interv
%             subplot(nr_interv+1,3,(i)*3+1)
%             topoplot(integrated(i,:),Info.Interface.EEG.chanlocs,'colormap',coolWarm);
%             colorbar;
%             a=get(gca,'xlim');
%             b=get(gca,'ylim');
%             text(a(2),b(1)-(b(2)-b(1))/10,'power (\muV^{2})')
%             
%             if strcmp(color_scale,'0-max')
%                 caxis([0 max(integrated(i,:))]);
%             else
%                 caxis([min(integrated(i,:)) max(integrated(i,:))]);
%             end
%             
%             text(-2,0.0525,[num2str(interval_Hz(i,1)),'-',num2str(interval_Hz(i,2)),'Hz'])
%         end
%         
%         %%               plot classiffication channels in 6 regions (to know which channels are in which regions)
%         
%         % figure(3)
%         % s{1}='r';
%         % s{2}='b';
%         % s{3}='k';
%         % s{4}='g';
%         % s{5}='y';
%         % s{6}='c';
%         % for i=1:6
%         %     plot(inty(ind{i}),intx(ind{i}),'.','color',s{i},'markersize',35)
%         %     hold on
%         % end
%         % hold off
%         % figure(1)
%         
%         %%           plot amplitudes of 129 channels in several frequency bands per channel
%         
%         for ii=1:nr_interv
%             subplot(nr_interv+1,3,(ii)*3+2)
%             if strcmp(color_scale,'0-max')
%                 nbt_plot_EEG_channels(integrated(ii,:),0 ,max(integrated(ii,:)),Info.Interface.EEG.chanlocs);
%             else
%                 nbt_plot_EEG_channels(integrated(ii,:),min(integrated(ii,:)) ,max(integrated(ii,:)),Info.Interface.EEG.chanlocs);
%             end
%             colorbar
%             a=get(gca,'xlim');
%             b=get(gca,'ylim');
%             text(a(2),b(1)-(b(2)-b(1))/10,'power (\muV^{2})')
%         end
%         
%         %%            for 6 regions, plot amplitudes averaged over all channels in region
%         
%         for ii=1:nr_interv
%             subplot(nr_interv+1,3,(ii)*3+3)
%             if strcmp(color_scale,'0-max')
%                 integrated_regions(ii,:)= nbt_plot_subregions(integrated(ii,:),1,0,max(integrated_regions(ii,:)));
%             else
%                 integrated_regions(ii,:)= nbt_plot_subregions(integrated(ii,:),1,min(integrated_regions(ii,:)),max(integrated_regions(ii,:)));
%             end
%             colorbar
%             axis equal
%             a=get(gca,'xlim');
%             b=get(gca,'ylim');
%             text(a(2),b(1)-(b(2)-b(1))/10,'power (\muV^{2})')
%         end
%         
%         
%         %%                            add Info to plot
%         
%         subplot(nr_interv+1,3,1)
%         text(0.5,0.5,'Interpolated topoplot','horizontalalignment','center')
%         axis off
%         
%         subplot(nr_interv+1,3,3);text(0.5,0.5,'Averaged in six subregions','horizontalalignment','center')
%         axis off
%         
%         subplot(nr_interv+1,3,2)
%         title(['Integrated amplitudes for experiment ',Info.file_name],'position',[0.5 1],'fontsize',11,'fontweight','bold','interpreter','none')
%         text(0.5,0.5,'Actual channels','horizontalalignment','center')
%         axis off
%         
%         
%         %%                   make PDF
%         
%         if save_pdf_plot
%             paperp=[0 0 30 25];
%             set(gcf,'paperunits','centimeters','papersize',[30 25],'paperposition',paperp)
%             print(gcf,'-dpdf',[Save_dir,'/',Info.file_name,' topoplot amplitudes.pdf'])
%             %             winopen([Info.file_name,' topoplot amplitudes.pdf'])
%         end
%         
%         %%
%         
%         %----------------Now do the same for the normalized amplitudes
%         
%         %--- get figure handle
%         
%         done=0;
%         figHandles = findobj('Type','figure');
%         for i=1:length(figHandles)
%             if strcmp(get(figHandles(i),'Tag'),'nz amplitudes')
%                 figure(figHandles(i))
%                 done=1;
%             end
%         end
%         if ~done
%             figure()
%             set(gcf,'Tag','nz amplitudes');
%         end
%         set(gcf,'numbertitle','off');
%         set(gcf,'name','NBT: Normalized amplitudes');
%         clf
%         
%         %---------plotting topoplots per frequency band
%         colormap(coolWarm);
%         for i=1:nr_interv
%             subplot(nr_interv+1,3,(i)*3+1)
%             topoplot(normalized(i,:),Info.Interface.EEG.chanlocs,'colormap',coolWarm);
%             colorbar
%             a=get(gca,'xlim');
%             b=get(gca,'ylim');
%             if strcmp(color_scale,'0-max')
%                 caxis([0 max(normalized(i,:))]);
%             else
%                 caxis([min(normalized(i,:)) max(normalized(i,:))]);
%             end
%             text(a(2),b(1)-(b(2)-b(1))/10,'power (%)')
%             text(-2,0.0525,[num2str(interval_Hz(i,1)),'-',num2str(interval_Hz(i,2)),'Hz'])
%         end
%         
%         %----plot amplitudes of 129 channels in several frequency bands per channel
%         
%         for ii=1:nr_interv
%             subplot(nr_interv+1,3,(ii)*3+2)
%             if strcmp(color_scale,'0-max')
%                 nbt_plot_EEG_channels(normalized(ii,:),0 ,max(normalized(ii,:)),Info.Interface.EEG.chanlocs);
%             else
%                 nbt_plot_EEG_channels(normalized(ii,:),min(normalized(ii,:)),max(normalized(ii,:)),Info.Interface.EEG.chanlocs);
%             end
%             colorbar
%             a=get(gca,'xlim');
%             b=get(gca,'ylim');
%             text(a(2),b(1)-(b(2)-b(1))/10,'power (%)')
%         end
%         
%         %%            for 6 regions, plot amplitudes averaged over all channels in region
%         
%         for ii=1:nr_interv
%             subplot(nr_interv+1,3,(ii)*3+3)
%             if strcmp(color_scale,'0-max')
%                 normalized_regions(ii,:)=nbt_plot_subregions(normalized(ii,:),1,0, max(normalized_regions(ii,:)));
%             else
%                 normalized_regions(ii,:)=nbt_plot_subregions(normalized(ii,:),1,min(normalized_regions(ii,:)), max(normalized_regions(ii,:)));
%             end
%             colorbar
%             axis equal
%             a=get(gca,'xlim');
%             b=get(gca,'ylim');
%             text(a(2),b(1)-(b(2)-b(1))/10,'power (%)')
%         end
%         
%         for i = 1:6
%             norm = sum(normalized_regions(:,i));
%             for j = 1:5
%                 normalized_regions(j,i) = normalized_regions(j,i) / norm;
%             end
%         end
%         
%         %%                         add Info to plot
%         
%         subplot(nr_interv+1,3,1)
%         text(0.5,0.5,'Interpolated topoplot','horizontalalignment','center')
%         axis off
%         
%         subplot(nr_interv+1,3,3);text(0.5,0.5,'Averaged in six subregions','horizontalalignment','center')
%         axis off
%         
%         subplot(nr_interv+1,3,2)
%         title(['Normalized amplitudes for experiment ',Info.file_name],'position',[0.5 1],'fontsize',11,'fontweight','bold','interpreter','none')
%         text(0.5,0.5,'Actual channels','horizontalalignment','center')
%         axis off
%         
%         %%                   make PDF
%         
%         if save_pdf_plot
%             paperp=[0 0 30 25];
%             set(gcf,'paperunits','centimeters','papersize',[30 25],'paperposition',paperp)
%             print(gcf,'-dpdf',[Save_dir,'/',Info.file_name,' topoplot normalized amplitudes.pdf'])
%             %             winopen([Info.file_name,' topoplot normalized amplitudes.pdf'])
%         end
%         
%     %end
end



%% update biomarker object 

amplitude_1_4_Hz =nbt_UpdateBiomarkerInfo(amplitude_1_4_Hz , Info);
amplitude_4_8_Hz =nbt_UpdateBiomarkerInfo(amplitude_4_8_Hz, Info);
amplitude_8_13_Hz =nbt_UpdateBiomarkerInfo(amplitude_8_13_Hz, Info);
amplitude_13_30_Hz =nbt_UpdateBiomarkerInfo(amplitude_13_30_Hz, Info);
amplitude_30_45_Hz =nbt_UpdateBiomarkerInfo(amplitude_30_45_Hz, Info);
amplitude_55_125_Hz = nbt_UpdateBiomarkerInfo(amplitude_55_125_Hz, Info);
amplitude_60_90_Hz = nbt_UpdateBiomarkerInfo(amplitude_60_90_Hz, Info);

amplitude_1_4_Hz_Normalized =nbt_UpdateBiomarkerInfo(amplitude_1_4_Hz_Normalized, Info);
amplitude_4_8_Hz_Normalized =nbt_UpdateBiomarkerInfo(amplitude_4_8_Hz_Normalized, Info);
amplitude_8_13_Hz_Normalized  =nbt_UpdateBiomarkerInfo(amplitude_8_13_Hz_Normalized, Info);
amplitude_13_30_Hz_Normalized=nbt_UpdateBiomarkerInfo(amplitude_13_30_Hz_Normalized, Info);
amplitude_30_45_Hz_Normalized=nbt_UpdateBiomarkerInfo(amplitude_30_45_Hz_Normalized, Info);
amplitude_55_125_Hz_Normalized=nbt_UpdateBiomarkerInfo(amplitude_55_125_Hz_Normalized, Info);
amplitude_60_90_Hz_Normalized=nbt_UpdateBiomarkerInfo(amplitude_60_90_Hz_Normalized, Info);
end

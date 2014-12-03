%------------------------------------------------------------------------------------
% Originally created by Simon-Shlomo Poil (2012), see NBT website for current email address
%------------------------------------------------------------------------------------
%
% ChangeLog - see version control log at NBT website for details.
%
% Copyright (C) 2012  Simon-Shlomo Poil 
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
% ---------------------------------------------------------------------------------------


function [ NewDataMatrix, BiomarkersToUse ] = nbt_RemoveFeatures( DataMatrix,outcome,Method,ChannelsToUse,UniqueBiomarkers)

NRofChannels = size(DataMatrix,2)/UniqueBiomarkers;


    

ChannelsToUse = ChannelsToUse(:);
switch Method
    case 'all'
        BiomarkersToUse = cell(1,1);
        if(isempty(ChannelsToUse))
            %dummy operation. Nothing to do.
            NewDataMatrix = DataMatrix;
            BiomarkersToUse{1,1} = 1:size(DataMatrix,2);
        else
            NewDataMatrix = nan(size(DataMatrix,1), UniqueBiomarkers);
            
            for Bid = 1:UniqueBiomarkers
                NewDataMatrix(:,Bid) = nanmedian(DataMatrix(:,(ChannelsToUse + NRofChannels*(Bid-1))),2);
                BiomarkersToUse{:,Bid} = (ChannelsToUse + NRofChannels*(Bid-1));
            end
        end
    case 'ttest2'
        % This largely ignores interaction effects
        if(isempty(ChannelsToUse)) %we prune accross all channles, without taking into account unique biomarkers.
            BiomarkersToUse = cell(1,1);
            p = NaN(size(DataMatrix,2),1);
            for iotta=1:size(DataMatrix,2)
                [h,p(iotta)] = ttest2(DataMatrix(outcome==1,iotta),DataMatrix(outcome==0,iotta),[],[],'unequal');
            end
            NewDataMatrix = DataMatrix(:, p<0.05);
            BiomarkersToUse{1,1} = find(p<0.05);
        else
            BiomarkersToUse = cell(1,UniqueBiomarkers);
            p = NaN(length(ChannelsToUse), 1);
            for Bid=1:UniqueBiomarkers
                for iotta=1:length(ChannelsToUse)
                    [h,p(iotta)] = ttest2(DataMatrix(outcome==1,(ChannelsToUse(iotta) + NRofChannels*(Bid-1))), DataMatrix(outcome==0,(ChannelsToUse(iotta) + NRofChannels*(Bid-1))));
                end
                if(~isempty(find(p<0.05)))
                    NewDataMatrix(:,Bid) = nanmedian(DataMatrix(:,(ChannelsToUse(p<0.05) + NRofChannels*(Bid-1))),2);
                    BiomarkersToUse{:,Bid} =  (ChannelsToUse(p<0.05) + NRofChannels*(Bid-1));
                else
                    NewDataMatrix(:,Bid) = nanmedian(DataMatrix(:,(ChannelsToUse(:) + NRofChannels*(Bid-1))),2);
                    BiomarkersToUse{:,Bid} =  (ChannelsToUse(:) + NRofChannels*(Bid-1));
                end
                if (sum(isnan(NewDataMatrix(:,Bid)))) %if nan we average across all channels
                    NewDataMatrix(:,Bid) = nanmedian(DataMatrix(:,(ChannelsToUse(:) + NRofChannels*(Bid-1))),2);
                    BiomarkersToUse{:,Bid} =  (ChannelsToUse(:) + NRofChannels*(Bid-1));
                end
            end
        end
    case 'ttest2-MCP'
        if(isempty(ChannelsToUse)) %we prune accross all channles, without taking into account unique biomarkers.
            BiomarkersToUse = cell(1,1);
            p = NaN(size(DataMatrix,2),1);
            for iotta=1:size(DataMatrix,2)
                [h,p(iotta)] = ttest2(DataMatrix(outcome==1,iotta),DataMatrix(outcome==0,iotta),[],[],'unequal');
            end
              p = nbt_MCcorrect(p,'bino');
            NewDataMatrix = DataMatrix(:, p);
            BiomarkersToUse{1,1} = p;
        else
            BiomarkersToUse = cell(1,UniqueBiomarkers);
            p = NaN(length(ChannelsToUse), 1);
            BBid = 0;
            for Bid=1:UniqueBiomarkers
                for iotta=1:length(ChannelsToUse)
                    [h,p(iotta)] = ttest2(DataMatrix(outcome==1,(ChannelsToUse(iotta) + NRofChannels*(Bid-1))), DataMatrix(outcome==0,(ChannelsToUse(iotta) + NRofChannels*(Bid-1))));
                end
                p = nbt_MCcorrect(p,'bino');
                if(~isempty(p))
                    BBid = BBid +1;
                    NewDataMatrix(:,BBid) = nanmedian(DataMatrix(:,(ChannelsToUse(p) + NRofChannels*(Bid-1))),2);
                    BiomarkersToUse{:,BBid} =  (ChannelsToUse(p) + NRofChannels*(Bid-1));
                        
                    if (sum(isnan(NewDataMatrix(:,BBid)))) %if nan we average across all channels
                        NewDataMatrix(:,BBid) = nanmedian(DataMatrix(:,(ChannelsToUse(:) + NRofChannels*(Bid-1))),2);
                        BiomarkersToUse{:,BBid} =  (ChannelsToUse(:) + NRofChannels*(Bid-1));
                    end
                end
            end
        end
    case 'ttest'
        %we assume the data is paired!
        if(isempty(ChannelsToUse)) %we prune accross all channles, without taking into account unique biomarkers.
            p = NaN(size(DataMatrix,2),1);
            for iotta=1:size(DataMatrix,2)
                [h,p(iotta)] = ttest(DataMatrix(outcome==1,iotta),DataMatrix(outcome==0,iotta),[],[],'unequal');
            end
            NewDataMatrix = DataMatrix(:, p<0.05);
            BiomarkersToUse{1,1} = find(p<0.05);
        else
            BiomarkersToUse = cell(1,UniqueBiomarkers);
            p = NaN(length(ChannelsToUse),1);
            for Bid=1:UniqueBiomarkers
                for iotta=1:length(ChannelsToUse)
                    [h,p(iotta)] = ttest(DataMatrix(outcome==1,(ChannelsToUse(iotta) + NRofChannels*(Bid-1))), DataMatrix(outcome==0,(ChannelsToUse(iotta) + NRofChannels*(Bid-1))));
                end
                if(~isempty(find(p<0.05)))
                    NewDataMatrix(:,Bid) = nanmedian(DataMatrix(:,(ChannelsToUse(p<0.05) + NRofChannels*(Bid-1))),2);
                    BiomarkersToUse{:,Bid} =  (ChannelsToUse(p<0.05) + NRofChannels*(Bid-1));
                else
                    NewDataMatrix(:,Bid) = nanmedian(DataMatrix(:,(ChannelsToUse(:) + NRofChannels*(Bid-1))),2);
                    BiomarkersToUse{:,Bid} =  (ChannelsToUse(:) + NRofChannels*(Bid-1));
                end
                if (sum(isnan(NewDataMatrix(:,Bid)))) %if nan we average across all channels
                    NewDataMatrix(:,Bid) = nanmedian(DataMatrix(:,(ChannelsToUse(:) + NRofChannels*(Bid-1))),2);
                    BiomarkersToUse{:,Bid} =  (ChannelsToUse(:) + NRofChannels*(Bid-1));
                end
            end
        end
    case 'SIS' % Sure Independence Screening
        % Fan Lv 2007 Sure Independence Screening for Ultra-High
        % Dimensional Feature Space
        
        BiomarkersToUse = cell(1,1);
        B_1=NaN(size(DataMatrix,2),1);
        for iotta=1:size(DataMatrix,2)
            B_1(iotta) = glmfit(DataMatrix(:,iotta),outcome,'binomial','link','logit','constant','off');
        end
        [B_1,IX]=sort(abs(B_1),'descend');
        biom_to_use_END=min(round(1/2*size(DataMatrix,1)),size(DataMatrix,1));%we select the lesser of the two:
        % either subset d_n=5.5*n^(2/3) for n initial samples OR the
        % initial predictor size p
        BiomarkersToUse{1,1}=IX(1:biom_to_use_END);
        BiomarkersToUse{1,1}=sort(BiomarkersToUse{1,1});
        NewDataMatrix = DataMatrix(:,BiomarkersToUse{1,1});
    case 'ISIS' % EXPERIMENTAL Iterative Sure Independence Screening
        %offers better performance in correlated biomarkers
        BiomarkersToUse = cell(1,1);
        B_1=NaN(size(DataMatrix,2),1);
        for iotta=1:size(DataMatrix,2)
            B_1(iotta) = glmfit(DataMatrix(:,iotta),outcome,'binomial','link','logit','constant','off');
        end
        [B_1,IX]=sort(abs(B_1),'descend');
        biom_to_use_END=min(round(5.5*size(DataMatrix,1)^(2/3)),size(DataMatrix,2));%we select the lesser of the two:
        biom_to_use_END=3;
        % either subset d_n=5.5*n^(2/3) for n initial samples OR the
        % initial predictor size p
        BiomarkersToUse{1,1}=IX(1:biom_to_use_END);
        BiomarkersToUse{1,1}=sort(BiomarkersToUse{1,1});
        unused_bioms=setdiff(1:size(DataMatrix,2),BiomarkersToUse{1,1});
        if biom_to_use_END~=size(DataMatrix,2)%only do this if there are predictors left
            [B,FitInfo] = lassoglm(DataMatrix(:,BiomarkersToUse{1,1}),outcome,'binomial',ones(biom_to_use_END,1),'Alpha',0.5,'DFmax',length(outcome),'CV',4,'NumLambda',50,'Standardize',1);
            lam=FitInfo.IndexMinDeviance;
            B0 = B(:,lam);
            cnst = FitInfo.Intercept(lam);
            B1 = [cnst;B0];
            preds = glmval(B1,DataMatrix(:,BiomarkersToUse{1,1}),'logit');
            residuals=abs((outcome - round(preds)));
            % Defining NEW biomarkers to add
            B_1NEW=NaN(size(DataMatrix(:,biom_to_use_END+1:end),2),1);
            for iotta=1:size(DataMatrix(:,biom_to_use_END:end),2)-1
                B_1NEW(iotta) = glmfit(DataMatrix(:,unused_bioms(iotta)),residuals,'binomial','link','logit','constant','off');
            end
            [B_1NEW,IXNEW]=sort(B_1NEW,'descend');
            biom_to_use_END_NEW=min(round(2.5*size(DataMatrix,1)^(2/3)),size(DataMatrix(:,biom_to_use_END+1:end),2));
            biom_to_use_END_NEW=3;
            bioms_to_add=unused_bioms(IXNEW(1:biom_to_use_END_NEW))';
            %             bioms_to_add=bioms_to_add(1:biom_to_use_END_NEW);
            BiomarkersToUse=[BiomarkersToUse{1,1};bioms_to_add];
            BiomarkersToUse{1,1}=sort(BiomarkersToUse{1,1});
            NewDataMatrix = DataMatrix(:,BiomarkersToUse{1,1});
        end
end
end


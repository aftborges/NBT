
classdef nbt_FeedbackSignal < nbt_Biomarker  % define here the name of the new object, here we choose nbt_Biomarker_template
    properties
        % add here the fields that are specific for your biomarker.
        %See the definition of nbt_Biomarker for fields that are allready there. For example:
        feedbackSignal
        ECRfeedbackSignal
        pBelowMedLine
        pAboveMedLine
        ECRpBelowMedLine
        ECRpAboveMedLine
    end
    methods
        % Now follows the definition of the function that makes a biomarker
        % of the type "nbt_Biomarker_template". The name of this function should alway be
        % the same as the name of the new biomarker object, in this example nbt_Biomarker_template
        % The inputs contain the information you want to add to the biomarker object :
        function BiomarkerObject = nbt_FeedbackSignal(feedbackSignal,ECRfeedbackSignal,pBelowMedLine,pAboveMedLine,ECRpBelowMedLine,ECRpAboveMedLine)
            
            % assign values for this biomarker object:
            BiomarkerObject.feedbackSignal = feedbackSignal;
            BiomarkerObject.ECRfeedbackSignal = ECRfeedbackSignal;
            BiomarkerObject.pBelowMedLine = pBelowMedLine;
            BiomarkerObject.pAboveMedLine = pAboveMedLine;
            BiomarkerObject.ECRpBelowMedLine = ECRpBelowMedLine;
            BiomarkerObject.ECRpAboveMedLine = ECRpAboveMedLine;
            
            %make list of biomarkers in this object:
            BiomarkerObject.Biomarkers ={'feedbackSignal','ECRfeedbackSignal','pBelowMedLine','pAboveMedLine','ECRpBelowMedLine','ECRpAboveMedLine'};
            
        end
        
        function plotThese(obj)
            figure
            plot(obj.feedbackSignal,'r');
            hold on
            plot(obj.ECRfeedbackSignal,'k');
        end
    end
    
end


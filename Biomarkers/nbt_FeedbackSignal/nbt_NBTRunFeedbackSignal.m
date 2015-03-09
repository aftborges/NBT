function nbt_runFeedbackSignal(Signal, SignalInfo, SaveDir)
feedbackSignal = nbt_doFeedbackSignal(Signal,SignalInfo,SaveDir);
nbt_SaveClearObject('feedbackSignal',SignalInfo,SaveDir);
end
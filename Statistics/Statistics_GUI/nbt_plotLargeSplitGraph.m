function mycallbackfcn(hObject, eventdata)
    hFig = get(hObject, 'parent');
    if strcmp(get(hFig,'SelectionType'),'open')
    % this was a double click, do something
        thisAx = hObject;
        figure;
        new_handle = copyobj(thisAx,gcf);
        set(gcf,'Position',[100 100 840 630]);
        set(new_handle,'Position',[0.1300    0.1100    0.7750    0.8150])
       
    end





end
function [nameRegex, bndry, evType, uType] = epoch_definitions()

nameRegex = {...
    '25-11_calibratie', ...
    '25-11_session', ...
    '26-11_calibratie', ...
    '26-11_session' ...
    };

bndry = {...
    % Boundaries for 25-11_calibratie
    [...
    165    220; ...
    232    287; ...
    321    434; ...
    454    569; ...
    707    867; ...
    874    995; ...
    1012   1126; ...
    1142   1257; ...
    1275   1322; ...
    1335   1388 ...
    ]; ...
    % Boundaries for 25-11_session
    [(0:500:2000)' (500:500:2500)']; ...
    % Boundaries for 26-11_calibratie
    [...
    446   502; ...
    510   563; ...
    664   770; ...
    799   914; ...
    937   1056; ...
    1070  1183; ...
    1205  1320; ...
    1336  1449; ...
    1466  1518; ...
    1527  1584; ...
    ]; ...
    % Boundaries for 26-11_session
    [(0:500:3000)' (500:500:3500)'] ...
    };

chunkTypes = cell(7,1);
for i = 1:numel(chunkTypes),
    chunkTypes{i} = ['chunk' num2str(i)];
end

evType = { ...
    % Boundaries for 25-11_calibratie
    { ...
    'eo'; ...
    'ec'; ...
    'excite'; ...
    'content'; ...
    'amuse'; ...
    'disgust'; ...
    'sadness'; ...
    'fear'; ...
    'eo'; ...
    'ec' ...
    }; ...
    % Boundaries for 25-11_session
    chunkTypes(1:5); ...
    % Boundaries for 26-11_calibratie
    { ...
    'eo'; ...
    'ec'; ...
    'excite'; ...
    'disgust'; ...
    'content'; ...
    'fear'; ...
    'amuse'; ...
    'sadness'; ...
    'eo'; ...
    'ec' ...
    }; ...
    % Boundaries for 26-11_session
    chunkTypes(1:7) ...
    };


uType = [];
for i = 1:numel(evType)
    uType = [uType;evType{i}];
end
uType = unique(uType);



end
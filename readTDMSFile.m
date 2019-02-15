function [d, p] = readTDMSFile(file)
% Read TDMS file and convert into HDF5 format
% Written CF
% Tested 8/4
% example call:
% file = 'F:\FlexibleWingProject\Programs\Software\TDMS_File_IO\test.tdms';
% d = readTDMSFile(file)

dat = TDMS_readTDMSFile(file);
for i = 1:length(dat.chanNames)
    for j = 1:length(dat.chanNames{i})
        var = dat.chanNames{i}{j};
        if ~strcmp(var,'Time')
            var = strrep(var,' ','_');
            ind = dat.chanIndices{i}(j);
            data = dat.data{ind};
            d.(var).Data = data;
            if ~isempty(dat.propValues{ind})
                props = dat.propValues{ind};
                t = double(props{3})*double(1:props{4});
            end
            d.(var).Axis = t;
        end
    end
end

%% Check if there are any params to extract from prop def
p.Stiffness = [];

for i =1:length(dat.propNames{1})
    types = {'Stiffness', 'Damping', 'Mass'};
    for j = 1:length(types)
        if ~isempty(strfind(dat.propNames{1}{i},types{j}))
            p.(types{j}) = dat.propValues{1}{i};
        end
    end
end
    
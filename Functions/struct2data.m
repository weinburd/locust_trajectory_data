function [data_final,data_neighbors] = struct2data(data_struct)
% Converts variables in the format of data_final and data_neighbors
% Into a single struct array

disp('Now converting your struct array into data_final')
tic

[Ntracks, Ntimes] = size(data_struct);
feats = vertcat(data_struct(:,1).features);
Nfeats = size(feats,2);

data_final = nan(Ntracks,Nfeats,Ntimes);
data_neighbors = cell(Ntracks,Ntimes);

for loc = 1:Ntracks
    for t = 1:Ntimes
        if ~isempty(data_struct(loc,t).features)
            data_final(loc,:,t) = data_struct(loc,t).features;
        end
        if ~isempty(data_struct(loc,t).neighbors)
            data_neighbors{loc,t} = data_struct(loc,t).neighbors;
        end
    end
end

% data = num2cell(data_struct);
% data = cellfun(@struct2cell, data, 'UniformOutput', false);
% 
% final_cell = cellfun(@(x) x{1}, data, 'UniformOutput', false);
% 
% data_final = cell2mat(final_cell);
% 
% data_neighbors = cellfun(@(x) x{2}, data, 'UniformOutput', false);
% 
% data_neighbors = cell(Ntracks,Ntimes);
% data_neighbors(:,:) = data(:,:){2};

fprintf(['That took %f seconds', newline],toc)

end

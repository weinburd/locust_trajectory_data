function data_new = data2struct(data_final,data_neighbors)
% Converts variables in the format of data_final and data_neighbors
% Into a single struct array

disp('Now converting your data into a struct array')
tic

[Ntracks, Nfeats, Ntimes] = size(data_final);

%indices that are not nan values
notnans = ~isnan( data_final(:,1,:) );
notnans = reshape( notnans, Ntracks, Ntimes);

%%% if we later want to divide the field features into subfields
% global  idx_x idx_y idx_flag...
%         idx_speed idx_theta...
%         idx_localSpeed idx_localStd idx_localMinMax...
%         idx_state
% 
% idx_x = 1; idx_y = 2; idx_flag = 3;
% idx_speed = 4; idx_theta = 5;
% idx_localSpeed = 6; idx_localStd = 7; idx_localMinMax = 8;
% idx_state = 9;
% 
% posn = makeFeatureCell(data_final,idx_x:idx_y,notnans);
% flag = makeFeatureCell(data_final,idx_flag,notnans);
% speed = makeFeatureCell(data_final,idx_flag,notnans);
% theta = makeFeatureCell(data_final,idx_flag,notnans);
% motionStats = makeFeatureCell(data_final,idx_flag,notnans);
% state = makeFeatureCell(data_final,idx_flag,notnans);

final_cell = makeFeatureCell(data_final, 1:Nfeats, notnans);

% notempty = find(~cellfun(@isempty,data_neighbors));
% final_neighbors = cell(Ntracks,Ntimes);
% final_neighbors(notempty) = data_neighbors(notempty);

data_new = struct( 'features', final_cell, 'neighbors', data_neighbors);

%%% Naive way with for loops %%%
% data_new = struct();
% 
% for loc = 1:Ntracks
%     for t = 1:Ntimes
%         
%         if ~isnan( data_final(loc,1,t) )
%             if sum( isnan(data_final(loc,1,t)) ) ~= 0
%                 error([ 'The x coord was NaN but at least one of the other features was not' newline...
%                         'Stopped at (locust, time) = (' loc ',' t ')']);
%                 % ^This never triggered, worth knowing!
%             end
%             data_new(loc,t).features = data_final(loc,:,t);
%         end
% 
%         if ~isempty( data_neighbors{loc,t} )
%             data_new(loc,t).neighbors = data_neighbors{loc,t};
%         end
%     end
% end

fprintf(['That took %f seconds', newline],toc)

%% Functions

function final_cell = makeFeatureCell(data_final,idx,notnans)
    data_this = data_final(:,idx,:);
    
    %reshape the data into a cell array
    data_cell = num2cell( data_this, 2);
    data_cell = reshape(data_cell, Ntracks, Ntimes);
    
    %ensure we don't encode a bunch of pointers to empty arrays
    final_cell = cell(Ntracks,Ntimes);
    final_cell(notnans) = data_cell(notnans);
end

end

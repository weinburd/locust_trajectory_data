function [data_Neighbors, neighbor_cloud] = dataNeighbors(data_final, fieldDims)
%%% Input:
% data_final = an Nlocust by 9 by Ntimestep matrix
%       each locust has 9 features  at each timestep
%       the features are:
%   [x y flag speed theta MagAveVelocity localStdSpeed min(fwdMax,bkwdMax) motionState]
% fieldDims = a 4 by 1 vector of the video's inscribed rectangular field dimentions (in physical units)
%       [left_side right_side bottom top]
%%% Outputs:
% data_Neighbors = an Nlocust by Ntimestep cell array
%                   each entry contains an array of locust indices that are neighbors of that locust at that timestep
% neighbor_cloud = an NtotalNeighbors by 2 array
%                  relative positions of all neighbors
disp('Now computing neighbors in your data');
tic
global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state

W = fieldDims(2)-fieldDims(1); % physical width of the video field
H = fieldDims(4)-fieldDims(3); % physical height of the video field

% extract position and orientation information
data = [data_final(:,idx_x:idx_y,:) data_final(:,idx_theta,:)];

[Nlocust, Nfeatures, Ntimesteps] = size(data);

neighbors = cell(Nlocust,Ntimesteps);
neighbor_cloud = cell(Ntimesteps,1);
%debugging
% global thetas
% thetas = [];
for t = 1:Ntimesteps
    d = data(:,:,t);
    [in_box, box_neighbors, posn_idx] = LocInBox(d);
    
    [neighbor_idx, neighbor_posn] = MakeNeighData(d, in_box, box_neighbors, posn_idx);
    
    neighbor_cloud{t} = neighbor_posn;
    neighbors(:,t) = neighbor_idx(:); 
end

neighbor_cloud = cell2mat(neighbor_cloud);
data_Neighbors = neighbors;
fprintf(['That took %f seconds', newline],toc)

%% Associated PiC functions

function [in_box, box_neighbors, posn_idx] = LocInBox(data)
    %%% Inputs:
    % data = N by 3 array. Rows = [xposn yposn headingdirection]
    %        (N = Nlocust)
    %%% Outputs:
    % in_box = Nx by Ny cell array. 
    %          Entries are vectors of locust indices contained in each box
    % box_neighbors = Nx by Ny cell array. 
    %                 Entries are vectors of locust indices contained in that box and all adjecent boxes.
    %                 Boxes on the edges are empty, since locusts in those boxes are counted as having no neighbors.
    % posn_idx = N by 1 cell array. Rows = { [xposition_index], [yposition_index] }

    N = size(data,1); %length(data); % number of agents

    dx = 7; %width of boxes (except ends)
    dy = 7; %height of boxes (except ends)

    % posn_idx = an N by 2 cell array
    % each row tells you the index of the box that locust is in
    % we use a centered grid, with extra space in outside boxes, 
    % for use with a video frame of size W x H

    Nx = round(W/dx); % = 14 %number of horizontal boxes
    Ny = round(H/dy); % = 8 %number of horizontal boxes

    shiftx = fieldDims(1) + W/2 - Nx/2*dx;
    shifty = fieldDims(3) + H/2 - Ny/2*dy; % for the 133 video this is negative
    % the effect is that we have boxes on the top/bottom with height < dy

    posn_idx = ceil( (data(:,1:2)-[shiftx,shifty])./[dx,dy] );

    posn_idx( posn_idx == 0 ) = 1; % extend small x,y edges
    posn_idx( posn_idx == Nx+1 ) = Nx; % extend right side box
    posn_idx( posn_idx(:,2)==Ny+1 , 2) = Ny; % extend large y box

    posn_idx = num2cell(posn_idx); % for indexing reasons later

    % an Nx by Ny cell array 
    % each entry is a list of locusts (indices) contained in that box
    in_box = cell(Nx,Ny);
    
    % Assign each locust's index to a box
%     for ii = 1:N
%         if ~isnan( data(ii,1) ) % if the locust IS in the frame
%             %debugging
%             test = posn_idx(ii,:);
%             if test{1} > Nx && test{2} > Ny
%                 disp('A locust is outside your particle in cell boxes, probably because of the way you scaled the width and height.')
%                 pause
%             end
% 
%             in_box{ posn_idx{ii,:} }(end+1) = ii;
%         end
%     end
    % slightly faster way to loop?
    idx_notNaN = ~isnan( data(:,1) );
    idx_notNaN = find(idx_notNaN);
    for ii = idx_notNaN'
        in_box{ posn_idx{ii,:} }(end+1) = ii;
    end

    % an Nx by Ny cell
    % each entry is a list of locusts (indices) contained in that box and adjacent boxes
    % NOT ANYMORE - However! We keep the edge entries empty, signifying that any locusts in an edge box have no neighbors (since they are not focal locusts).
    box_neighbors = cell(Nx, Ny);

    %%% TO DO: 
    %%% DONE - Modify idx_list so that we can iterate over all cells, i.e. x = 1:Nx
    %%% DONE - Then include distance from edge in the min(dx,dy) when in complex form
    %%% This allows us to extract more data, still unbiased by the edge.
    for x = 1:Nx
        for y = 1:Ny
            % create an index mask: 25 boxes with box (x,y) at the center
            idx_list1 = { x-1, y-1; x, y-1; x+1, y-1;
                         x-1, y;   x, y;   x+1, y;
                         x-1, y+1; x, y+1; x+1,y+1  };
            idx_list2 = { x-2, y-2; x-1, y-2; x, y-2; x+1, y-2; x+2, y-2;
                          x-2, y-1;                             x+2, y-1;
                          x-2, y;                               x+2, y;
                          x-2, y+1;                             x+2, y+1;
                          x-2, y+2; x-1, y+2; x, y+2; x+1, y+2; x+2, y+2  };
            idx_list = [idx_list1; idx_list2];
            % deal with edge cases
            idx_rmx = zeros(length(idx_list),1); idx_rmy = zeros(length(idx_list),1);
            if x <= 2
                idx_rmx = (cell2mat(idx_list(:,1))<=0); %left edge
            end
            if x >= Nx-1
               idx_rmx = (cell2mat(idx_list(:,1))>=Nx+1); %right edge
            end
            if y <= 2
               idx_rmy = (cell2mat(idx_list(:,2))<=0); %top edge
            end
            if y >= Ny-1
               idx_rmy = (cell2mat(idx_list(:,2))>=Ny+1); %bottom edge
            end
            idx_rm = idx_rmx | idx_rmy; %this OR statement ensures no double removing
            idx_list = idx_list(~idx_rm,:);
            
            theseNeighbors = [];
            for idx = 1:length(idx_list)
                theseNeighbors = [theseNeighbors in_box{ idx_list{idx,:} }];
            end
            % This appears to hurt performance for fewer than 500 locusts,
            %   also for 593 locust as in 133_30sec_1120_1150
            % OR preallocate as a cell array for speed
%             theseNeighbors = cell(1,length(idx_list));
%             for idx = 1:length(idx_list)
%                 theseNeighbors{idx} = in_box{ idx_list{idx,:} }; % when preallocating
%             end
%             theseNeighbors = cell2mat(theseNeighbors); % convert to an array, when preallocating
            

            box_neighbors{x,y} = theseNeighbors;
        end
    end

end

function [neighbor_idx, neighbor_posn] = MakeNeighData(data, in_box, box_neighbors, posn_idx)
    %%% Inputs:
    % data = N by 3 array. Rows = [xposn yposn headingdirection]
    % in_box = Nx by Ny cell array. 
    %          Entries are vectors of locust indices contained in each box
    % box_neighbors = Nx by Ny cell array. 
    %                 Entries are vectors of locust indices contained in that box and all adjecent boxes.
    %                 Boxes on the edges are empty, since locusts in those boxes are counted as having no neighbors.
    % posn_idx = N by 1 cell array. Rows = { [xposition_index], [yposition_index] }
    % countedge = scalar. For debugging
    % 
    %%% Outputs:
    % neighbor_idx = N by 1 cell where each row is for a locust and the entry
    %               on that row is an array of indices for that locust's
    %               neighbors
    % neighbor_posn = unknown by 2 array with relative positions of all locusts
    %               in the frame
    
    N = length(data); % number of agents

    dx = 7; %width of boxes (except ends)
    dy = 7; %height of boxes (except top and bottom)
    Nx = size( in_box , 1); %14;
    Ny = size( in_box, 2); %8;

    neighbor_posn = cell(N,1); % preallocate as a cell array for speed
    neighbor_idx = cell(N,1);
    % Set to NaN if locust is not in the frame
    %nans = isnan(data(:,1));
    %neighbor_idx( nans,1 )= {NaN};

    for x = 1:Nx
        for y = 1:Ny
            for loc = in_box{x,y}
%                 edgeflag = posn_idx{loc,1}==1 | posn_idx{loc,1}==Nx |...
%                            posn_idx{loc,2}==1 | posn_idx{loc,2}==Ny ;
%                 edgeflag = false; %this was when we weren't allowing any focals on the edges.
%                 if edgeflag % if the locust is around the edge
%                     %neighbor_idx{loc} = NaN; %maybe it's easier to store with nothing in all those cells?
%                 else
                poss_neigh = box_neighbors{x,y};

                rel_posn = data(poss_neigh,1:2) - data(loc,1:2);

                dist = sqrt(rel_posn(:,1).^2+rel_posn(:,2).^2);
                dedge = edgedist(data(loc,1:2)); %distance from focal locust to the edge
                
                %we consider only neighbors that are:
                %   a) not the same as the focal locust
                %   b) within 2 box widths from the focal locust
                %   c) within the distance from the edge to the focal locust (removes side to side edge biases)
                close_neigh = 0 < dist & dist < min([2*dx,2*dy,dedge]);
                these_neighs = poss_neigh( close_neigh );
                
                %if the locust has no neighbors, we can save some space by
                %leaving that element of the cell array as Null, rather
                %than explicitly assigning it the empty array [] as a value.
                if ~isempty(these_neighs)
                    neighbor_idx{loc} = these_neighs;
                end
                
                % convert to complex numbers
                % throw out the neighbor positions that are not close enough
                %rel_posn = rel_posn( close_neigh, : );
                rel_posn(~close_neigh,:) = [];
                
                % y-component gets a negative sign due to ImageJ pixel coordinate
                crel_posn = rel_posn(:,1)-1i*rel_posn(:,2); %complexify
                % rotate by focal locust's orientation data(loc,3)
                % and by a factor of pi/2 so that focal locust is facing upwards
                crel_posn = abs(crel_posn)...
                            .*exp(1i*angle(crel_posn) - 1i*data(loc,3) + 1i*pi/2);
                
                %debugging
                %if ~isempty(poss_neigh( close_neigh ))
                %    thetas = [thetas data(loc,3)];
                %end
                
                % convert back to (x,y) positions
                rel_posn = [real(crel_posn) imag(crel_posn)];

                neighbor_posn{loc} = rel_posn;
%                 end
            end
        end
    end

    % convert to an array
    neighbor_posn = cell2mat(neighbor_posn);
    
    %%% A different (better?) way to index boxes?
    % for box_idx = 1:numel(in_box) % could exclude edges to speed up
    %     for locust = in_box{box_idx}
    %         edgeflag = posn_idx{locust,1}==1 | posn_idx{locust,1}==Nx |...
    %                    posn_idx{locust,2}==1 | posn_idx{locust,2}==Ny ;
    %         if edgeflag % if the locust is around the edge
    %             neighbor_idx{locust} = NaN;
    %             edgecount = edgecount + 1;
    %         else
    %         neigh_idx = box_neighbors{box_idx};
    %         
    %         rel_posn = data(neigh_idx,1:2)-data(locust,1:2);
    %         
    %         % convert to complex numbers
    %         crel_posn = rel_posn(:,1)+1i*rel_posn(:,2);
    %         % consider only neighbors within one box width
    %         neigh_idx = neigh_idx( abs(crel_posn) < min(dx,dy) );
    %         
    % %         % rotate by focal locust's orientation, again considering only neighbors within one box
    % %         crel_posn = crel_posn(  0 < abs(crel_posn) &...
    % %                                 abs(crel_posn) < min(dx,dy) )*exp(-1i*data(locust,3));
    % %         % convert back to (x,y) positions
    % %         rel_posn = [real(crel_posn) imag(crel_posn)];
    % %         
    %         neighbor_posn{locust} = 1; %rel_posn;
    %         
    %         neigh_idx = neigh_idx( neigh_idx ~= locust );
    %         neighbor_idx{locust} = neigh_idx;
    %         end
    %     end
    % end

end

function distToEdge = edgedist(posn)
    xFoc = posn(1);
    yFoc = posn(2);
    left = xFoc-fieldDims(1); right = fieldDims(2)-xFoc;
    bottom = yFoc-fieldDims(3); top = fieldDims(4)-yFoc; % is this right? pixels vs. y vals, etc...
    distToEdge = min([left, right, bottom, top]);
end

end
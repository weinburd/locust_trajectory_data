function [thisData, varargout] = neighPositions(data, neighbors, varargin)
%%% Input:
% data = an Nlocust by 4 by Ntimestep matrix
%           for each locust has four features at each timestep
%           the features are [x_posn y_posn orientation motionState]
% neighbors = an Nlocust by Ntimestep cell array
%           each entry contains an array of locust indices that are neighbors of that locust at that timestep
% varargin{1} = inDebug = struct with four fields = 
%           inDebug.count = (total neighbors)
%           inDebug.countnan = (focal not in frame?), 
%           inDebug.countempty= (focal w/no neighs),
%           inDebug.countnoneighs = (focal w/no neighs?)}
% varargin{2} = reshapeData = 1 by 2 cell array = {reshapeOption, bodyEllipse(assumed shape)}
%
%%% Output:
% thisData = nNeighs by 3 matrix = [xpos ypos focalState]
% varargout{1} = outDebug = 1 by 4 cell array = adjusted values for 

global numNeighbors numFocal

nOptInputs = numel(varargin);
if nOptInputs > 1
    inDebug = varargin{1};
    reshapeData = varargin{2};
elseif nOptInputs > 0
    inDebug = varargin{1};
    reshapeData = {'none', [0 0] };
else
    inDebug.count = 0;
    inDebug.countnan = 0;
    inDebug.countempty = 0;
    inDebug.countnoneighs = 0;
    reshapeData = {'none', [0 0] };
end

reshapeOpt = reshapeData{1};
% when reshapeOpt = 'rescale', rescales all x,y coords.
% when reshapeOpt = 'subtract', subtracts an ellipse from all x,y coords.
bodyEllipse = reshapeData{2};

Ntimesteps = size(neighbors,2);
thisData = cell(Ntimesteps,1);

    for t = 1:Ntimesteps
        % find the nonempty entries for this frame, i.e. the indices of locusts that have neighbors
        locusts = find( ~cellfun( @isempty, neighbors(:,t) ) );
        nowNeighbors = cell(size(locusts,1),1);
        for locust = locusts' % locusts must be a ROW vector
            neigh_idx = neighbors{locust,t};
            if isnan(neigh_idx) % focal locust is either not in the frame or is in an edge box
                % This never happens anymore because I took out NaNs for locusts out of the frame or on the edge.
                %allData = cat(1,allData,[NaN, NaN]); %unnesseccary
                inDebug.countnoneighs = inDebug.countnoneighs + 1;
            else
                numFocal = numFocal + 1;
                rel_posn = data(neigh_idx,1:2,t)-data(locust,1:2,t);
        
                % convert to complex numbers
                % y-component gets a negative sign due to ImageJ pixel coordinate
                crel_posn = rel_posn(:,1)-1i*rel_posn(:,2); 
                
                
                % rotate by focal locust's orientation data(locust,3,t)
                % and by a factor of pi/2 so that focal locust is facing up
                crel_posn = abs(crel_posn)...
                            .*exp(1i*angle(crel_posn) - 1i*data(locust,3,t) + 1i*pi/2);
                
                %debugging
                %global THETAS
                %THETAS = [THETAS data(locust,3,t)];

                %neighbor orientation relative to focal locust orientation
                %rel_theta = wrap2Pi(data(neigh_idx,3,t)-data(locust,3,t));
                
                % reshape complex positions
                haxis = bodyEllipse(1)/2;
                vaxis = bodyEllipse(2)/2;
                if strcmp(reshapeOpt,'rescale')
                    % rescales y coords to make the ellipse into a circle with the same width
                    crel_posn = real(crel_posn) + 1i*haxis/vaxis.*imag(crel_posn);
                elseif strcmp(reshapeOpt, 'subtract')
                    % this finds the distance from the origin to the point
                    % on the ellipse at an angle = angle(crel_posn)
                    % NOTE that angle is NOT equal to the ellipse parameter.
                    hdistsqr = 1./(1/haxis^2+tan(angle(crel_posn)).^2/vaxis^2);
                    ell_dist = sqrt(hdistsqr + tan(angle(crel_posn)).^2.*hdistsqr);
                    sub_dist = max( abs(crel_posn)-ell_dist, 0);
                    % add a line to remove any folks at zero...?
                    crel_posn = (sub_dist)./abs(crel_posn).*crel_posn;
                elseif strcmp(reshapeOpt,'none')
                    % do nothing
                end

                if isa(numNeighbors,'double')
                    % consider only the prescribed number of neighbors
                    [~,I] = mink(abs(crel_posn), numNeighbors);
                    crel_posn = crel_posn(I);
                    last_idx = min(numNeighbors,length(crel_posn));
                elseif isa(numNeighbors,'char')
                    % consider all neighbors
                    last_idx = length(crel_posn);
                end

                % convert back to (x,y) coords
                rel_posn = [real(crel_posn) imag(crel_posn)];
                
                %%% debugging... %%%
                if isnan(rel_posn(1:last_idx,:)) %(data(locust,3,t)) %using data catches some places where the heading is undefined but there are no neighbors
                % focal locust was in the frame, but heading = NaN because it either just entered or just left the frame
                    inDebug.countnan = inDebug.countnan + size(rel_posn(1:last_idx,:),1);
                elseif sum(sum(isnan(rel_posn(1:last_idx,:))))>0 %if there are ANY NaN values
                    error('There is a NaN hiding in rel_posn, but it is not all NaNs')
                end
                if isempty(rel_posn(1:last_idx,:)) %(neigh_idx)
                % focal locust was in the frame and not in an edge box, but has no neighbor within min(dx,dy)
                    inDebug.countempty = inDebug.countempty + 1;
                end
                inDebug.count = inDebug.count + max(1,size(rel_posn(1:last_idx,:),1));
                % should have: count == countnan+countempty+size(rmmissing(allData),1)
                %        also: count == countempty+size(allData,1)
                
                theseNeighbors = rel_posn(1:last_idx,:);
            end
            nNeighs = size(theseNeighbors,1);
            nowNeighbors{locusts==locust} = [ theseNeighbors data(locust,4,t)*ones(nNeighs,1) ];
        end
        thisData{t} = cell2mat(nowNeighbors);
    end

thisData = cell2mat(thisData);

outDebug = inDebug; % adjusted in place
varargout{1} = outDebug;

end
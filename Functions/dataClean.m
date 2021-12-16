function data_clean = dataClean(data_rough, scale, trans, varargin)
%%% Input:
% data_rough = a matrix in the format that lines up with the input for Data_Cleaning.m
%               each row is a track
%               each track has [timestep xposition yposition ]
%               spatial units are pixels
% scale = scalar = scaling factor in pix/cm
% tform = 3 by 3 matrix = transformation from board corners in pixels to 
%           actual board dimensions in cm
% optional input 1 is a flag to import the data including motion state data
%%% Output:
% data_clean = an Nlocusts by 3 by Ntimesteps matrix
%               each locust has [x y flag] at each timestep
%                   OR
%               each locust has [x y motion_state flag]
%               the spatial units are cm
%               the flag values are:
%                   0 = absent data because locust is not in video's 
%                       field of view (not yet entered or already has exited)
%                   1 = data from when a locust is detected
%                   2 = missing data to be interpolated (track had a missing spot)
disp('Now cleaning your data');
tic

nInputs = numel(varargin);
if nInputs > 0
    manualflag = varargin{1};
else
    manualflag = false;
end

if manualflag
    numFeatures = 4;
else
    numFeatures = 3;
end

Minit = data_rough; 

if isempty(trans)
    % Use 1/1.94 for Jacob's manual data
    %xscale = scale; % Jacob had: (1/18.4); 
    %yscale = scale; % Jacob had: (1/18.4);

    % convert to centimeters
    Minit(:,2:numFeatures:end)= Minit(:,2:numFeatures:end) / scale;
    Minit(:,3:numFeatures:end)= Minit(:,3:numFeatures:end) / scale;
else
    % Apply the projective transformation trans, given as input by a 3x3 matrix
    xcoords = Minit(:,2:numFeatures:end);
    ycoords = Minit(:,3:numFeatures:end);
    
    projPoints = [xcoords; ycoords; ones(size(xcoords))];
    projPoints = reshape(projPoints, size(xcoords,1), []); %convert to projective coordinates
    bigTrans = kron(eye(size(xcoords,2)), sparse(trans));
    projNewPoints = projPoints*bigTrans; %apply transformation
    
    %convert back to Euclidian coords
    newxCoords = projNewPoints(:,1:3:end)./projNewPoints(:,3:3:end);
    newyCoords = projNewPoints(:,2:3:end)./projNewPoints(:,3:3:end);
    
    Minit(:,2:numFeatures:end) = newxCoords;
    Minit(:,3:numFeatures:end) = newyCoords;
end


%%% My version of obtaining the *reshaped* 3D Mnext matrix
frames = Minit(:,1:numFeatures:end);
Ntimes = max(frames(:))+1; %'omitnan' is default, frames start from 0 so we add +1
Ntracks = size(Minit, 1);

Mnext = nan(Ntracks, numFeatures, Ntimes);
for col = 1:numFeatures:size(Minit, 2)
    theseFrames = Minit(:,col)+1;
    idx = find( ~isnan(theseFrames) );
    theseFrames( isnan(theseFrames) ) = []; %remove NaNs
    
    if length( theseFrames ) ~= length(idx)
        error('length problem')
    end
    
    
    for i = 1:length(idx)
        Mnext( idx(i),1:numFeatures-1,theseFrames(i) ) = Minit( idx(i), col+1:col+numFeatures-1);
    end
    % can be vectorized?
%     newM = reshape( Minit(idx, col+1:col+numFeatures-1), length(idx), numFeatures-2, [] );
%     Mnext(idx, 1:numFeatures-1, theseFrames) = newM;
end

%%% Now we apply flags:
%                   0 = absent data because locust is not in video's 
%                       field of view (not yet entered or already has exited)
%                   1 = data from when a locust is detected
%                   2 = missing data to be interpolated (track had a missing spot)
for loc = 1:size(Mnext, 1)
    trackx = Mnext(loc,1,:);
    trackx = trackx(:);
    
    tracked = ~isnan(trackx);
    
    %set flag = 1 for a locust that was successfully tracked
    Mnext(loc,numFeatures,tracked) = 1; 
    
    first = find( tracked, 1, 'first');
    %set flag = 0 for locusts who have not entered the video field
    Mnext(loc,numFeatures,1:first-1) = 0; 
    last = find( tracked, 1, 'last');
    %set flag = 0 for locusts who have exited the video field
    Mnext(loc,numFeatures,last+1:end) = 0; 
    
    missed = logical([ zeros(first-1,1) ; ~tracked(first:last); zeros(Ntimes-(last),1) ]);
    %set flag = 2 for locusts who have entered the video field, but then go untracked before leaving
    Mnext(loc,numFeatures, missed ) = 2;  
end

%%% Now we use linear interpolation to fill in position coordinates
%%% wherever we have a missing detection (flag = 2).

Mfinal3D = Mnext;

flags = reshape( Mfinal3D(:,numFeatures,:), Ntracks, Ntimes );
[locusts, times] = find( flags == 2);
pairs = [locusts times];
pairs = sortrows(pairs, 1);
Npairs = size( pairs, 1);

clear idx
for idx = 1:Npairs
    locust = pairs(idx,1);
    time = pairs(idx,2);
    
    j = 1; %this will be the size of the gap, i.e. number of missing detections
    % while the locust is the same and the next time is the original + gapsize
    while locust == pairs(min(idx+j,Npairs),1) && time+j == pairs(min(idx+j,Npairs),2)
        j = j+1; %increase gapsize
    end
    
    pos1 = Mfinal3D(locust,1:2,time-1);
    pos2 = Mfinal3D(locust,1:2,time+j);
    Mfinal3D(locust,1:2,time) = linInterp(pos1,pos2,j,time);
end
    
data_clean = Mfinal3D;
fprintf(['That took %f seconds', newline],toc)


%%% Other Functions %%%

    function intPos = linInterp(pos1,pos2,j,time)
        intPos = pos1 + (time - (time-1))/(j+1) * (pos2-pos1);
    end

end
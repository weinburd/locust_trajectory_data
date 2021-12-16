function data_rough = xml2mat(file,varargin)
%%% Inputs:
% file = an .xml file of tracks only
%       obtained by dropdown option: "Save as .xml" --> Execute in TrackMate
%       or by the scripted tracking load_track.py or load_track.m
% optional input 1 = true produces a plot
% optional input 2 is a flag to treat the .xml as manually tracked data
nInputs = numel(varargin);
if nInputs > 1
    plotflag = varargin{1};
    manualflag = varargin{2};
elseif nInputs > 0
    plotflag = varargin{1};
    manualflag = false;
else
    plotflag = false;
    manualflag = false;
end
%%% Output:
% mat = a matrix in the format that lines up with the input for Data_Cleaning.m
%       each row is a track
%       each track has [timestep xposition yposition (manual color, if enabled)]
%       spatial units are pixels

if manualflag
    fprintf(['Now importing spots from %s', newline],file)
    tic
    [spot_table, spot_ID_map] = trackmateSpots(file);
    fprintf(['That took %f seconds', newline],toc)

    fprintf(['Now importing edges %s', newline],file)
    edge_features = {'EDGE_TIME','MANUAL_COLOR'};
    tic
    edge_map = trackmateEdges_JWB(file,edge_features);
    fprintf(['That took %f seconds', newline],toc)

    track_names = edge_map.keys;
    n_tracks = numel(track_names);

    tracks = cell(n_tracks,1);

    for s = 1:n_tracks
        track_name = track_names{s};
        edge_table = edge_map( track_name );
        edge_table = sortrows( edge_table, 'EDGE_TIME'); % sort the table by time

        spot_IDs = unique( [edge_table.SPOT_SOURCE_ID...
                            edge_table.SPOT_TARGET_ID] );
        spot_rows = cell2mat( values(spot_ID_map, num2cell(spot_IDs)) );
        spot_features = spot_table( spot_rows,{'POSITION_T',...
                                               'POSITION_X',...
                                               'POSITION_Y'}  );
        spot_features = sortrows(spot_features, 'POSITION_T');
        
        track = [];
        track(:,1:3) = table2array(  spot_features  );
        track(1,4) = NaN; % the first timestep doesn't have a well-defined motion state
        track(2:end,4) = table2array( edge_table( :, 'MANUAL_COLOR') );

        tracks{s} = track;
    end
    % Each track is now an N_det x 4 matrix, one line per detection (so track was
    % detected in N_det frames)
    % On each row, the spot data is:
    %       [timestep   x   y  manual_color]
    numFeatures = size(tracks{1},2);
elseif ~manualflag || isempty(varargin{2})
    
    clipZ = true; % cut off the z value, making each track a N_det x 3 matrix
    % Each track will be a N_det x 4 matrix, one line per detection (so track was
    % detected in N_det frames)
    % On each row, the spot data is:
    %       [timestep   x   y ]

    fprintf(['Now importing %s', newline],file)
    tic
    tracks = importTrackMateTracks(file,clipZ);
    fprintf(['That took %f seconds', newline],toc)
    numFeatures = size(tracks{1},2);
    
    n_tracks = numel( tracks );
end

fprintf(['Found %d tracks in the file.', newline], n_tracks)

if plotflag
    % plot from the cell array
    figure() % make a new figure
    for ii = 1:n_tracks
        hold on
        plot(tracks{ii}(:,2), tracks{ii}(:,3))
    end
end

%takes tracks as an import in cell array format and produces a matrix in
%the format that lines up with Data_Cleaning.m
fprintf(['Now reformatting tracks from %s', newline],file)

reformatMat=cell2mat(tracks);
TT = max(reformatMat(:,1));
%set up a matrix of NaN values in which to load the data
Import = NaN(length(tracks),(TT+1)*numFeatures);

for locust = 1:length(tracks)
    focalLoc = transpose(cell2mat(tracks(locust)));
    focalLoc = focalLoc(:);
    focalLoc = transpose(focalLoc);
    Import(locust,1:length(focalLoc)) = focalLoc;
end

data_rough = Import;

end
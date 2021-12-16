% This is a script to test the full data pipeline
% .xml --> data_final
% OR
% video --> data_final
% using functions

addpath('Functions',...
            'Functions/packages/CircStat',...
            'Functions/packages/TrackMateScripts',...
        'Data/examples', 'Data/manual')

%%% options %%%
plotFlag = false;
saveData = false;
computeAccuracy = false;

newData = true;
    newManual = false;
    newAutomatic = true; % imports from existing .xml data files if newTracking = false
        %suboptions
        newTracking = false; % NOTE: This requires an old version of TrackMate, see ReadMe
            %subsuboptions
            saveVid = true; %save a copy of the processed video
            %saveXML = true; %save a copy of the .XML data file, necessary to progress...
            %always true else, nothing to read new data from

% data processing options
smoothFlag = true;
interpThetaFlag = true;            
            
file_idx = 1;
recording_examples = initData(file_idx,'examples');

vid_Landsberg = '133_10sec_710_720';
vid_Sharma = '133_10sec_1910_1920';

vidList = {vid_Landsberg, vid_Sharma};

%% Processing Data %%
if newData
    
    %%% Spatial Scale Factors %%%
    % for the 133 videos: 1260 pi / 68 cm = 18.53 pi/cm, see image scale.tif
    scale133 = 18.34113771; %18.53; %pix/cm 18.34113771;
    
    scale_Sharma = 1/10*scale133; % for Tarush's manual data
    scale_Landsberg = 1.94 * 18.34113771/18.4; % for Jacob's manual data, because of how he obtained his downsampled video
                       % 1.94 = the scaling Jacob used for his manual data
                       % 18.4 = the scaling Jacob used for his automatic data
                       % 18.34... = the scaling we're using for our automatic tracking
    fieldDims = [0 1824 0 1026] / scale133;
    %fieldDims_manual = [182 102] / scale_manual;
    % W = fieldDims(1) = 98.4350, for Video 33 (originally was 1824/18.4 = 99.0761)
    % H = fieldDims(2) = 55.3697, for Video 33 (originally was 1026/18.4 = 55.7065)

    %%% initialize data variable %%%
    data = cell( 4, 4 );
    
    %%%%%%%%%% MANUAL DATA %%%%%%%%%%
    if newManual
    disp('Now handling your MANUAL data...');
    manualFlag = true;
    
    for vid_idx = 1:numel(vidList)
        idx = 2*vid_idx-1; %maps 1 -> 1 and 2 -> 3
        
        vid = vidList{vid_idx};
        data{idx+1,1} = ['manual_' vid '_motion'];
        
        if vid_idx == 1
            scale_manual = scale_Landsberg;
        elseif vid_idx == 2
            scale_manual = scale_Sharma;
        end
        fieldDims_manual = [0 182 0 102] / scale_manual; %fieldDims;%
        
        file_manual = [data{idx+1,1} '.xml'];%['manual_' vid '_motion.xml'];

        data_rough_manual = xml2mat(file_manual, plotFlag, manualFlag);
        
        data_clean_manual = dataClean(data_rough_manual, scale_manual, [], manualFlag);

        model = [];
        data_final_manual = dataForm(data_clean_manual,model, smoothFlag, interpThetaFlag, plotFlag, manualFlag);
        data{idx+1,2} = data_final_manual; % store the final data

        [data_neighbors_manual, neighbor_cloud_manual] = dataNeighbors(data_final_manual, fieldDims_manual);
        data{idx+1,3} = data_neighbors_manual; % store the neighbor data
        data{idx+1,4} = neighbor_cloud_manual; % store the neighbor cloud

        recording_examples(file_idx).data{6+idx+1,2} = data2struct(data_final_manual,data_neighbors_manual);
        recording_examples(file_idx).data{6+idx+1,3} = neighbor_cloud_manual;
    end
    end
    
    if newAutomatic
    %%%%%%%%%% AUTOMATIC DATA %%%%%%%%%%
    disp('Now handling your AUTOMATIC data...');
    manualFlag = false;
    if newTracking
        ImageJ %maybe with a check and possibly an optional input "false"

        %tracking parameters in cm
        radius = 0.523413550009270; % = 9.6 pix/(scale133 pix/cm)
        %maximal distance allowed when linking two spots to continue a track
        searchRadius = 4.307257338617954; % = 79.0 pix/(scale133 pix/cm)
        %maximal distance allowed when linking two spots to create a new track
        linkDistance = 2.017323057327396; % = 37.0 pix/(18.53 pix/cm)
        %minimal distanced necessary to keep a track
        %we set this to zero for the very short 10-second clips
        displacementFilter = 0; %2*2; % = 2-10 cm; approximately 1-5 body lengths
        %tracking parameters in pixels
        rad = radius * scale133;
        frameGap = 1;
        searchRad = searchRadius * scale133;
        linkDist = linkDistance * scale133;
        dispFilter = displacementFilter * scale133;
        trackParams = [rad, frameGap, searchRad, linkDist, dispFilter];

        disp('Now tracking locusts in the .AVI video');
        
        for vid_idx = 1:numel(vidList)
            vid = vidList{vid_idx};
            idx = 2*vid_idx-1;
            
            saveXML = true;
            file_scripted = track(vid, trackParams, saveVid, saveXML);
            data{idx,1} = file_scripted;
            
            % requires additional dependency package...
            % not necessary for small examples
            % very helpful for large batch tracking
            %jheapcl % clear out some of the Java Heap Memory
        end
        
        ij.IJ.run("Quit","");
    end
    
    %%% Train an SVM for motion state assignment %%%
    if ~exist('model','var')
        % use 'gaussian' or an empty function call for the most accurate SVM
        % use 'linear' for faster model fitting
        model = fitSVM('linear'); % fit the SVM using Jacob's data
    end
    
    for vid_idx = 1:numel(vidList)
    idx = 2*vid_idx-1;
    
    if vid_idx == 1
        scale_manual = scale_Landsberg;
    elseif vid_idx == 2
        scale_manual = scale_Sharma;
    end
    
    file_scripted = ['tracks_' vidList{vid_idx} '.xml']; %data_accuracy{idx,1};
    data{idx,1} = file_scripted;
    data_rough = xml2mat(file_scripted,plotFlag);

    data_clean = dataClean(data_rough, scale133, []);

    data_final = dataForm(data_clean,model,smoothFlag, interpThetaFlag, plotFlag);
    data{idx,2} = data_final; %store the final data

    [data_neighbors, neighbor_cloud] = dataNeighbors(data_final, fieldDims);
    data{idx,3} = data_neighbors;
    data{idx,4} = neighbor_cloud;
    
    recording_examples(file_idx).data{6+idx,2} = data2struct(data_final,data_neighbors);
    recording_examples(file_idx).data{6+idx,3} = neighbor_cloud;
    end
    end
    
else
    load('data_recording_examples.mat')
end


%% Save Data & Compute Accuracy%%

if saveData
    datafile = 'Data/examples/data_recording_examples.mat';
    save(datafile,'recording_examples');
end

data_final_L = struct2data(recording_examples(1).data{7,2});
data_final_L_manual = struct2data(recording_examples(1).data{8,2});
data_final_S = struct2data(recording_examples(1).data{9,2});
data_final_S_manual = struct2data(recording_examples(1).data{10,2});

%%% Compute accuracy metrics %%%
if computeAccuracy
unassignDist = 0.75; % in cm
unassignDist_link = sqrt(2)*unassignDist;

% need to offset data because they don't start at the exact same frame
% offset indices
auto_Lands = 1:249;
man_Lands = 2:250;

auto_Sharm = 1:247;
man_Sharm = 4:250;
auto_Sharm = 1:247; % By looking at the videos, I thought it should be 1:248 and 3:250
man_Sharm = 4:250;  % But this gives best values SSC ~ 0.80 and LSC ~ 0.80, when unassignDist = 0.5

disp('Now computing accuracy...');

SSC_L = spotAccuracy(data_final_L(:,:,auto_Lands), data_final_L_manual(:,:,man_Lands), unassignDist);
[LSC_L,FP,FN] = linkAccuracy(data_final_L(:,:,auto_Lands), data_final_L_manual(:,:,man_Lands), unassignDist_link);

SSC_S = spotAccuracy(data_final_S(:,:,auto_Sharm), data_final_S_manual(:,:,man_Sharm), unassignDist);
[LSC_S,FP,FN] = linkAccuracy(data_final_S(:,:,auto_Sharm), data_final_S_manual(:,:,man_Sharm), unassignDist_link);
    
fprintf(['Landsberg: SSC = %f | LSC = %f \n'...
         '   Sharma: SSC = %f | LSC = %f \n'],...
        SSC_L, LSC_L, SSC_S, LSC_S)
end

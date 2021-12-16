function [trackfile, new_minQuality]  = track(vid, trackParams, saveVid, saveXML, varargin)
%%% A MATLAB script for running our automated tracking procedure in ImageJ
%%% dependencies:
%       1) Fiji distribution of ImageJ
%       2) ImageJ-MATLAB project (add update site in Fiji) https://imagej.net/scripting/matlab
%       
%%% Input:
%       vidList = cell array containing strings of the names of video files to track
%       saveVid = a flag that chooses whether or not to save a .AVI of the processed video
%       saveXML = a flag that chooses whether or not to save the data to a new .XML file
%%% Output:
%   The function has no output, but executes the following operations...
%       *proccesses the video vid (ImageJ2 aka Fiji), optionally saves a copy
%       *automatically tracks locusts in vid (TrackMate
%       *shows the user a histogram plot of the Quality, and takes user input for a minimal threshold value
%       *repeats automatic tracking and applies the manually-input threshold value
%       *saves the data to a .XML file that can be imported to Matlab with xml2mat.m

% put this outside the function
% Initialize an instance of ImageJ2 (Fiji)
% ImageJ
% ij.IJ.run("Quit","");

%----------------------------------
% Import Fiji and TrackMate classes
%----------------------------------

import java.io.File
import java.util.HashMap
import ij.*

% import fiji.plugin.trackmate.TrackMate
% import fiji.plugin.trackmate.Model
% import fiji.plugin.trackmate.Settings
% import fiji.plugin.trackmate.SelectionModel
% import fiji.plugin.trackmate.Logger
% import fiji.plugin.trackmate.features.FeatureFilter
% import fiji.plugin.trackmate.detection.LogDetectorFactory
% import fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory
% import fiji.plugin.trackmate.tracking.LAPUtils
% import fiji.plugin.trackmate.gui.displaysettings.DisplaySettingsIO
% import fiji.plugin.trackmate.gui.displaysettings.DisplaySettings
% import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer

import fiji.plugin.trackmate.*
import fiji.plugin.trackmate.detection.*
import fiji.plugin.trackmate.features.*
import fiji.plugin.trackmate.features.track.*
import fiji.plugin.trackmate.tracking.*
import fiji.plugin.trackmate.visualization.hyperstack.*
import fiji.plugin.trackmate.providers.*
import fiji.plugin.trackmate.action.*
import fiji.plugin.trackmate.io.*

import fiji.stacks.Hyperstack_rearranger.*

% hacking things together to run without a loop but leaving the structure
% so I can possibly put the loop back in the future.
vidList = {vid};
vid = 1;
% for vid = 1:numel(vidList)
    %---------------------------------------
    % Point at a new video file
    %---------------------------------------
    IJ.run("Close All")
    
    proc_vid = strcat('preprocessed_',vidList{vid},'.avi');
    % This path depends on the machine you are working on...
    % Since these paths are passed to ImageJ, they must be defined as
    % absolute paths (rather than relative).
    % Therefore you must input a specific path for your instance code.
    if ispc
        % PC path
        savePath = 'C:\\Users\\jazzy\\Documents\\GitRepos\\locust_trajectory_data\\Data\\examples\\';
    elseif ismac
        % Mac
        savePath = '/Users/weinburd/Documents/GitRepos/locust_trajectory_data/Data/examples/';
    else
        error("You are not working on a Mac or a PC, so we didn't know where to find your video files.")
    end
	
    vidFile_in = [savePath, proc_vid];
    vidFile_out = [savePath, 'processed_', vidList{vid},'.avi' ];

    if isfile(vidFile_out)
        imp = IJ.openImage(vidFile_out);%open the new image
        % could also use command ImagePlus()?
        imp.show() %show/select image
    else
        imp = IJ.openImage(vidFile_in); %open the new image
        % could also use command ImagePlus()?
        imp.show() %show/select image
        
        %---------------------------------------
        % Video Processing in ImageJ (macro language)
        %---------------------------------------
        
        m8bit = 'run("8-bit"); '; 	%Fiji/Trackmate thinks our video is still in color. 
        %This has something to do with the codec we're using to save the video in an nv12 format.
        mInvert = 'run("Invert", "stack"); ';
        mBackSub = ['run("Z Project...", "projection=Median"); '...
                    'imageCalculator("Subtract create stack", '...
                    '"' proc_vid '", '...
                    '"MED_' proc_vid '");'];
        mBrightCon = ['run("Brightness/Contrast..."); '...
                          'setMinAndMax(0, 100); '...
                          'run("Apply LUT", "stack");'];
        mBlur = 'run("Gaussian Blur...", "sigma=4 stack");';
        %mReorder = 'run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Slices (z)] frames=[Frames (t)]");';
        %macro = strcat(m8bit, mInvert, mBackSub, mBrightCon, mBlur, mReorder);
        % Adding this to the macro creates a Java Heap Memory error
        
        macro = strcat(m8bit, mInvert, mBackSub, mBrightCon, mBlur);
        
        IJ.runMacro(macro); %
        
        % Save the processed image
        if saveVid
            %vidFile_out = [savePath, vidList{vid},'.avi' ];
            mSaveAVI = strcat('run("AVI... ", "compression=None frame=25 save=', vidFile_out, '");');
            IJ.runMacro(mSaveAVI);
            disp(['Saved a video file to' newline vidFile_out])
        end
        
        imp = WindowManager.getCurrentImage();
        
    end % if isfile(vidFile_out)

    % switch Z slices and T frames
    imp = reorderHyperstack(imp,'CTZ',false,true); 
    % false => leave old image open (elminates a user dialogue box to save)
    % true => open new image

	%imp.show() %unnecessary
	disp('The image link will be')
	disp(imp)
    
    % requires additional dependency package...
    % not necessary for small examples
    % very helpful for large batch tracking
    %jheapcl % clear out some of the Java Heap Memory
    
    %----------------------------
    % Create the model object now
    %----------------------------

    % Some of the parameters we configure below need to have
    % a reference to the model at creation. So we create an
    % empty model now.
    model = Model();

    % Send all messages to ImageJ log window.
    model.setLogger(Logger.IJ_LOGGER)
    logger = model.getLogger();
    
    %---------------------------
    % Unpack tracking parameters
    %---------------------------
    rad = trackParams(1);
    frameGap = int32(trackParams(2)); %2
    %int32() ensure Jave reads it as an integer. See -- https://www.mathworks.com/help/matlab/matlab_external/data-type-conversions.html
    searchRad = trackParams(3); %70.0
    linkDist = trackParams(4); %10.0
    dispFilter = trackParams(5); % 4-10 cm
    %------------------------
    % Prepare settings object
    %------------------------
    
    %settings = Settings(imp); %new version of TrackMates
    settings = Settings();
    settings.setFrom(imp) %old version of TrackMate
    
    % Configure detector - We use a java map
    settings.detectorFactory = LogDetectorFactory();
    map = HashMap();
    map.put('DO_SUBPIXEL_LOCALIZATION', true);
    map.put('RADIUS', rad);
    map.put('TARGET_CHANNEL', 1);
    map.put('THRESHOLD', 0.1);
    map.put('DO_MEDIAN_FILTERING', false);
    settings.detectorSettings = map;
    
    %settings.detectorSettings.put('RADIUS',rad);
        
    % Configure tracker
    settings.trackerFactory  = fiji.plugin.trackmate.tracking.kalman.KalmanTrackerFactory();
    map2 = HashMap();
    map2.put('MAX_FRAME_GAP',frameGap); %2
    map2.put('KALMAN_SEARCH_RADIUS',searchRad); %70.0
    map2.put('LINKING_MAX_DISTANCE',linkDist); %10.0
    settings.trackerSettings = map2;

    %-------------------
	% Instantiate plugin
	%-------------------
	
	trackmate = TrackMate(model, settings);
	       
	%--------
	% Process
	%--------
	    
    ok = trackmate.checkInput();
    if ~ok
        display(trackmate.getErrorMessage())
    end

    ok = trackmate.process();
    if ~ok
        display(trackmate.getErrorMessage())
    end
	
    %----------------------------------------
	% Take a User-Input minumum Quality
    % and rerun tracking with a new filter
	%----------------------------------------
    
    % 1) Get spot quality from Model
	spots = model.getSpots();
    qualities = spots.collectValues('QUALITY', false);
    % 2) plot histogram using
    figure(31)
    histogram(qualities)
    % 3) ask for user input value
    nOptInput = numel(varargin);
    if nOptInput == 0
        new_minQuality = input('Choose a new minimum quality value \n');
    else
        new_minQuality = 0.95;
    end
    % 4) apply a filter with that value for Quality
    %       (An alternative would be to set the threshold in the detector,
    %       which could optimize things for bigger videos if needed.)
    % Configure spot filters - Classical filter on quality
	% In the user interface this value gets chosen when you look at the histogram of Qualities.
% 	settings.clearSpotFilters() %to remove the old filter
%     filter1 = FeatureFilter('QUALITY', new_minQuality, true);
% 	settings.addSpotFilter(filter1)
    %do this instead
    settings.detectorSettings.put('THRESHOLD',new_minQuality);
    % Configure track filters.
    % The displacement feature is provided by the TrackDurationAnalyzer.
    settings.addTrackAnalyzer(TrackDurationAnalyzer())
    filter2 = FeatureFilter('TRACK_DISPLACEMENT', dispFilter, true);
    settings.addTrackFilter(filter2)
    
    % filter out tracks with short duration
    %filter3 = FeatureFilter('TRACK_DURATION', 10, true);
    %settings.addTrackFilter(filter3)
    
    % 5) do "Process" above again...
	% Instantiate Plugin
    trackmate = TrackMate(model, settings);
	% Process  
    ok = trackmate.checkInput();
    if ~ok
        display(trackmate.getErrorMessage())
    end

    ok = trackmate.process();
    if ~ok
        display(trackmate.getErrorMessage())
    end
    
    %----------------
	% Display results
	%----------------
	
	selectionModel = SelectionModel(model);
	displayer =  HyperStackDisplayer(model, selectionModel, imp);
	displayer.render()
	displayer.refresh()
	    
	% Echo results with the logger we set at start:
	logger.log('-----------New settings and New Model-----------')
	logger.log(string(model))
    logger.log(string(settings))
    
    disp('Inspect the tracking and press a key to continue!')
    pause
    
	%----------------
	% Export the resulting model to a .xml file
	%----------------
    
	trackfile=['tracks_', vidList{vid}, '.xml'];
    if saveXML
        
        file_out = File([savePath, trackfile]);
        % This is NOT the full Trackmate .xml file. 
        % It is just the data .xml file obtained in the GUI by "Save Tracks as .xml" -> "Execute".
		ExportTracksToXML.export( model, settings, file_out)
    end
% end

%%% 

IJ.run("Close All")

% put this outside the function
% Quit the current instance of ImageJ2
% ij.IJ.run("Quit","");

end
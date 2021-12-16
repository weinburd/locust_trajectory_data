function SSC = spotAccuracy(testauto, testman, costUnmatched)
%%% Inputs:
%%%     testauto = matrix of data from automatic tracking
%%%     testman = matrix of data from manual tracking
%%%     **Must have same number of timesteps = size(input, 3)
%%%     **Timesteps must line up perfectly
%%%     costUnmatched = the cost of not matching a spot
%%%         **This is the cost associated with leaving a spot in 'auto' NOT matched to
%%%         one in 'manual'. Since our costs are just Euclidian distance, this
%%%         basically says the following:
%%%         For a spot X in 'auto', if there is no spot in 'manual' that is within a
%%%         distance of costUnassigned centimeters from X, then leave X unassigned.

%%% Outputs:
%%%     JSC = Jacob Similarity Coefficient =  TP/(TP + FP + FN);

% isolate the position data
testauto = testauto(:,1:2,:);
testman = testman(:,1:2,:);

%%Here we assemble a Cost matrix of distances from the auto to the manual
%%points. For example, Cost(i,j) is the Euclidean distance from point i to point j
%%The Cost matrix has the number of rows equal to the number of tracks in
%%the automatic tracks and the number of columns equal to the number of
%%tracks in the manual track. Its depth is the number of frames.
Cost = zeros(length(testauto),length(testman),size(testauto, 3));
for j = 1:size(testauto, 3)
    for k = 1:length(testauto)
        Cost(k,:,j) = vecnorm( [testauto(k,1,j)-transpose(testman(:,1,j));...
                                testauto(k,2,j)-transpose(testman(:,2,j))]...
                                ,2,1)';
    end
end

%%Now we set all NaN values to a super high value in the cost matrix, so
%%they will not be matched with anything. So long as they are higher than
%%the unassignment cost, I think this will be fine. This is just because
%%matchpairs can't handle NaN values. Match is a matrix where each row
%%shows the indices of paired points in automatic and manual track. So if a
%%row is [3,4], then the point in the 3rd row of testauto corresponds to
%%the point in the 4th row of testman.
Cost(isnan(Cost))=1000;

%unassignment = .5; % this is a distance (in cm)... beyond which points won't be matched?
unassignment = costUnmatched;

%for m = 0:.01:1.25
    %unassignment = m;
    Match = NaN(length(testman),2,size(testauto,3));
    for j = 1:size(testauto, 3)
        singleMatch = matchpairs(Cost(:,:,j),unassignment);
        Match(1:size(singleMatch,1),:,j) = singleMatch;
    end

    %%Now I'll make a matrix to find the number of points present in both the
    %%manual and automatic tracks at each time step. Row 1 will be the number
    %%in the automatic, Row 2 will be the number in the manual, Row 3 will
    %%be the number of paired spots, Row 4 will be the number of false
    %%positives (Row 1 minus Row 3), andRow 5 will be the number of false
    %%negatives (Row 2 minus Row 3).
    PointsPresent = zeros(5,1,size(testauto,3));
    for j = 1:size(testauto, 3)
        PointsPresent(1,1,j) = sum(~isnan(testauto(:,1,j)));
        PointsPresent(2,1,j) = sum(~isnan(testman(:,1,j)));
        PointsPresent(3,1,j) = sum(~isnan(Match(:,1,j)));
        PointsPresent(4,1,j) = PointsPresent(1,1,j)-PointsPresent(3,1,j);
        PointsPresent(5,1,j) = PointsPresent(2,1,j)-PointsPresent(3,1,j);
    end
    
    %%JS will be our version of the Jaccard Similarity Coefficient. It will be
    %%equal to the true positives, i.e., the total number of matched points
    %%within the unassignment distance (TP), divided by the sum of the true
    %%positives plus the number of false positives (FP) plus the number of
    %%false negatives (FN).
    TP = sum(PointsPresent(3,1,:),3);
    FP = sum(PointsPresent(4,1,:),3);
    FN = sum(PointsPresent(5,1,:),3);
    SSC = TP/(TP + FP + FN);
end
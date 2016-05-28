%% Ricardo Sanchez Castillo
%Id Number 4225015
%User psxrs4@nottingham.ac.uk
%% Overview
% There are a numerous algorithms for binocular stereo, specifically for 
% image matching. However, they can be classified into two main categories: 
% Sparse matching algorithms where feature matching is based on the 
% strongest features such as edges or corners, and although this benefits 
% the performance of an algorithm the output is only a sparse disparity 
% data. On the other hand, are the dense matching algorithms mainly based 
% on block matching along the rectified epipolar lines of the two stereo 
% images; however this approach is subjected to many errors leading to 
% noisy disparity maps. This algorithm is based on Zhang and Shan, 
% A Progressive Scheme for Stereo Matching, using both feature and 
% template matching. 
%% Parameters
warning('off','signal:findpeaks:largeMinPeakHeight')
testImage = 'teddy';
firstImageFile = '-im2.png';
secondImageFile = '-im6.png';
SURFMethod = 'MetricThreshold';
SURFMethodValue = 600;
maxRatioFeatureMatching = 0.4;
matchFeaturesMethod = 'MaxRatio';
matchedFeaturesThreshold = 0.25;
radius = 2;
maxIterations = 5000;
disparity = 1;
variance_max = 1.0;
variance_min = 0.3;
variance_control = 30;
correlation_window_sizex = 11;
correlation_window_sizey = 9;
minPeakHeight = 0.70;

UNKNOW = 0;
MATCHED = 1;
NOMATCH = 2;

%% Step 1 Feature extraction.
% SURF features will be obtained from both images using image pyramids, 
% but only the 1500 strongest features points will be used for matching.

firstImage = rgb2gray(imread(strcat(testImage, firstImageFile)));
firstImagePyramid1 = impyramid(firstImage, 'reduce');
firstImagePyramid2 = impyramid(firstImagePyramid1, 'reduce');

secondImage = rgb2gray(imread(strcat(testImage, secondImageFile)));
secondImagePyramid1 = impyramid(secondImage, 'reduce');
secondImagePyramid2 = impyramid(secondImagePyramid1, 'reduce');

imagePoint1L1 = detectSURFFeatures(firstImagePyramid1, SURFMethod, SURFMethodValue);
imagePoint2L1 = detectSURFFeatures(secondImagePyramid1,SURFMethod,SURFMethodValue);
imageFeatures1L1 = extractFeatures(firstImagePyramid1, imagePoint1L1);
imageFeatures2L1 = extractFeatures(secondImagePyramid1, imagePoint2L1);

imagePoint1L2 = detectSURFFeatures(firstImagePyramid2, SURFMethod, SURFMethodValue);
imagePoint2L2 = detectSURFFeatures(secondImagePyramid2, SURFMethod, SURFMethodValue);
imageFeatures1L2 = extractFeatures(firstImagePyramid2, imagePoint1L2);
imageFeatures2L2 = extractFeatures(secondImagePyramid2, imagePoint2L2);

%% Step 2 Feature matching
% The features will be matched using the matchFeatures function provided 
% by Matlab which has been seen to function with enough accuracy.

pairsFirstLevel = matchFeatures(imageFeatures1L1, imageFeatures2L1, matchFeaturesMethod, maxRatioFeatureMatching);
pairsSecondLevel = matchFeatures(imageFeatures1L2, imageFeatures2L2, matchFeaturesMethod, maxRatioFeatureMatching);

% Remove those points where the Y level are not equal in both matches due 
% to the epipolar restriction
pairsMatched = [];
for r = 1 : size(pairsFirstLevel, 1)
    first = imagePoint1L1(pairsFirstLevel(r, 1)).Location;
    second = imagePoint2L1(pairsFirstLevel(r, 2)).Location;
    if abs(first(:, 2) - second(:, 2)) < matchedFeaturesThreshold
        pairsMatched = vertcat(pairsMatched, pairsFirstLevel(r, :)); %#ok<AGROW>
    end
end
pairsFirstLevel = pairsMatched;
pairsMatched = [];
for r = 1 : size(pairsSecondLevel, 1)
    first = imagePoint1L2(pairsSecondLevel(r, 1)).Location;
    second = imagePoint2L2(pairsSecondLevel(r, 2)).Location;
    if abs(first(:, 2) - second(:, 2)) < matchedFeaturesThreshold
        pairsMatched = vertcat(pairsMatched, pairsSecondLevel(r, :)); %#ok<AGROW>
    end
end
pairsSecondLevel = pairsMatched;

matchedPointFirstLevel1 = imagePoint1L1(pairsFirstLevel(:, 1));
matchedPointFirstLevel2 = imagePoint2L1(pairsFirstLevel(:, 2));
matchedPointSecondLevel1 = imagePoint1L2(pairsSecondLevel(:, 1));
matchedPointSecondLevel2 = imagePoint2L2(pairsSecondLevel(:, 2));
leftMatchedPoints = matchedPointFirstLevel1.Location * 2;
leftMatchedPoints = vertcat(leftMatchedPoints, matchedPointSecondLevel1.Location * 4);
rightMatchedPoints = matchedPointFirstLevel2.Location * 2;
rightMatchedPoints = vertcat(rightMatchedPoints, matchedPointSecondLevel2.Location * 4);
leftMatchedPoints = round(leftMatchedPoints);
rightMatchedPoints = round(rightMatchedPoints);
leftMatchedPoints = horzcat(leftMatchedPoints, repmat(variance_min, [size(leftMatchedPoints, 1) 1]));

%% Step 3 Pixel labelling
% All pixels except those which correspond to SURF features will be 
% labelled as UNKNOWN.
pixels = zeros(size(firstImage));
for r = 1: size(leftMatchedPoints, 1)
    pixels(leftMatchedPoints(r, 2), leftMatchedPoints(r, 1)) = MATCHED;
end
%% Step 4 Compute correlation on pixels around the seeds.
% Zhang and Shan algorithm propose to consider first those pixels in 
% highly textured regions followed by those with a large uncertainty.
% However, for this algorithm is proposed start with pixels around those
% that have been labelled as MATCHED, therefore for the first iteration,
% pixels around the SURF features will be considered first and ratio will
% be incremented each iteration. The correlation technique to be used will
% be the same as the original algorithm used by Zhang and Shan.
hasPixelsMatched = 1;
statistics = [];
while hasPixelsMatched == 1
    %For each matched point
    hasPixelsMatched = 0;
    for r = 1 : size(leftMatchedPoints, 1)
        x = leftMatchedPoints(r, 1);
        y = leftMatchedPoints(r, 2);
        x1 = rightMatchedPoints(r, 1);
        y1 = rightMatchedPoints(r, 2);
        %Take the pixels in the neightbourhood with ratio r
        for s = -radius : radius
            for c = -radius : radius
                if x + c < 1 || x + c > size(pixels, 2) || y + s < 1 || y + s > size(pixels, 1)
                    continue
                end
                %We only take UNKNOW pixels
                if pixels(y + s,x + c) ~= UNKNOW || (s == 0 && x == 0)
                    continue
                end
                %Candidate pixels
                f_k = sqrt((disparity * (s ^ 2)) / (1 - ((disparity ^ 2) / 4)));
                x_min = round(c + x1 - f_k);
                if x_min < 0
                    x_min = 0;
                end
                x_max = round(c + x1 + f_k);
                if x_max > size(firstImage, 2)
                    x_max = size(firstImage, 2);
                end
                xs = (x_min : x_max)';
                candidatePixels = horzcat(xs, repmat(y + s, [size(xs, 1) 1]));
                if size(candidatePixels, 1) == 0
                    continue
                end
                correlation = zeros(size(candidatePixels, 1), 1);
                                
                %Compute correlation based on candidate pixels
                %Template
                temp_x_min = x + c - correlation_window_sizex;
                if temp_x_min <= 0
                    temp_x_min = 1;
                end
                temp_y_min = y + s - correlation_window_sizey;
                if temp_y_min <= 0
                    temp_y_min = 1;
                end
                temp_x_max = x + c + correlation_window_sizex;
                if temp_x_max > size(firstImage, 2)
                    temp_x_max = size(firstImage, 2);
                end
                temp_y_max = y + s + correlation_window_sizey;
                if temp_y_max > size(firstImage, 1)
                    temp_y_max = size(firstImage, 1);
                end
                    
                template = firstImage(temp_y_min : temp_y_max, temp_x_min : temp_x_max);
                %Compute correlation coefficient for each candidate pixels
                for i = 1 : size(candidatePixels, 1)
                    temp_x_min = candidatePixels(i, 1) - correlation_window_sizex;
                    if temp_x_min <= 0
                        temp_x_min = 1;
                    end
                    temp_y_min = candidatePixels(i,2) - correlation_window_sizey;
                    if temp_y_min <= 0
                        temp_y_min = 1;
                    end
                    temp_x_max = candidatePixels(i, 1) + correlation_window_sizex;
                    if temp_x_max > size(secondImage, 2)
                        temp_x_max = size(secondImage, 2);
                    end
                    temp_y_max = candidatePixels(i,2) + correlation_window_sizey;
                    if temp_y_max > size(secondImage, 1)
                        temp_y_max = size(secondImage, 1);
                    end
                   corr_win = secondImage(temp_y_min : temp_y_max, temp_x_min : temp_x_max); 
                   if size(template, 1) < size(corr_win, 1) || size(template, 2) < size(corr_win, 2)
                       corr_win = corr_win(1:size(template, 1), 1:size(template, 2));
                   end
                   if size(corr_win, 1) < size(template, 1) || size(corr_win, 2) < size(template, 2)
                       template = template(1:size(corr_win, 1), 1:size(corr_win, 2));
                   end
                   correlation(i, 1) = corr2(template, corr_win);
                end               
                %Normalize correlation to be between 0 and 1
                correlation = correlation + 1;
                correlation = correlation / 2;
                %Detect peaks
                if size(candidatePixels, 1) == 1 && correlation(1,1) > minPeakHeight
                    peaks = correlation(1, 1);
                elseif size(candidatePixels, 1) == 2 && max(correlation) > minPeakHeight
                    peaks = max(correlation);
                elseif size(candidatePixels, 1) > 3
                    peaks = findpeaks(correlation, 'MinPeakHeight', minPeakHeight);
                else
                    peaks = [];
                end
                
                %If there is one peak label as MATCHED
                if size(peaks, 1) == 1 && peaks(1,1) < 1
                    element = find(correlation == peaks);
                    e = candidatePixels(element, :);
                    pixels(y + s,x + c) = MATCHED;
                    leftMatchedPoints = vertcat(leftMatchedPoints, [x + c, y + s, 0]); %#ok<AGROW>
                    rightMatchedPoints = vertcat(rightMatchedPoints, e(1, :)); %#ok<AGROW>
                    hasPixelsMatched = 1;
                end
                
                %If there is any peak label as NOMATCH
                if size(peaks, 1) == 0
                    pixels(y + s,x + c) = NOMATCH;
                    hasPixelsMatched = 1;
                end
            end
        end
    end
    %Stores statistics
    statistics = vertcat(statistics, [size(leftMatchedPoints, 1), size(find(pixels==0), 1), size(find(pixels==2), 1)]); %#ok<AGROW>
    if mod(size(leftMatchedPoints, 1), 1000) == 0
        size(leftMatchedPoints)
    end
end

%% Step 4 Create disparity map
% Shows the results, calculating the disparity in each pixel and showing as
% at image
disparity = leftMatchedPoints(:, 1) - rightMatchedPoints(:, 1);
pixels1 = pixels;
for r = 1 : size(leftMatchedPoints, 1)
    pixels1(leftMatchedPoints(r, 2), leftMatchedPoints(r)) = disparity(r);
end
I = mat2gray(pixels1);
figure
subplot(2, 2, 1), imshow(firstImage)
subplot(2, 2, 2), imshow(I)
subplot(2, 2, 3), imshow(secondImage)
subplot(2, 2, 4), imshowpair(firstImage, I)
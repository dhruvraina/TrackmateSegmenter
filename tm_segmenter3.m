% TrackMate segmenter:
% Author: d.raina
% last edit: 11.April.17

%PseudoCode:

%Read in trackmate 'matlab' file. Use each point as a "seed" point to
%perform nuclear segmentation. Read out median intensities. Currently this
%works by using matlab's 'MultiThresh' applied locally. KMSeg is also an
%option to implement!

close all; clear all
%ToDo: File i/o ??  ::c_dir should be passed through the datloader-11/04/17
addpath('/Users/draina/Documents/MATLAB/Dhruv''s Stuff/MPI_Progs/ExchangeTools/DrosteEffect-BrewerMap-a77e675')
addpath('/Users/draina/Documents/MATLAB/Dhruv''s Stuff/MPI_Progs/TrackmateTools')

%% Inputs:
% c_dir = '/Users/draina/Desktop/SiR_Test/';
% file = [c_dir 'RegularOutput.xlsx'];

c_dir = '/Users/draina/Desktop/Draq5_Test/';
file = [c_dir 'InvertedTrack_matlab.xml'];
[trackList trackMeta] = importTrackMateTracks(file);

% rawTrackImport = readtable(file);

mainImg = tiffread2([c_dir '5mMDRAQ_Washout_200_inv1']);
measurementImg = tiffread2([c_dir '5mMDRAQ_Washout_200_inv1']);
scalingFactor = 1;  %:: Scaling factor is derived from image Metadata. TBD: automatically extract this number.


%% For Inverted Tracks:
invertmode = 'y';
framestart = 1;
frameend = length(mainImg);
if invertmode =='y'
    
    %Make a backup of the original variable for debugging
    templist = trackList;
    
    %Convert the frame numbers into reverse order
    gg = cellfun(@(x) frameend-x(:,1), trackList, 'UniformOutput', 0);
    
    %Replace reversed frame numbers in original trackList
    for cc = 1:length(trackList)
        trackList{cc,1}(:,1)= gg{cc,1}(:,1);
    end
    
    %Re-sort accordind to new frame numbers:
    trackList = cellfun(@(x) sortrows(x,1), trackList, 'UniformOutput', 0);
    
end



%% Framewise Segmentation and Datextraction:
for ctr1 = 1:length(mainImg)
%for ctr1 = 140:230
  % ctr1 = 201;
    
    %Find all tracks which are active in this frame:
    c_img          = mainImg(ctr1).data;
    c_measureImg   = measurementImg(ctr1).data;
    c_imgTrax1     = cellfun(@(x)find(x==ctr1),trackList(:,1), 'UniformOutput', false);
    c_imgTrax      = cellfun(@(x)~isempty(x), c_imgTrax1, 'UniformOutput', false);
    [c_imgTrax, q] = find(cell2mat(c_imgTrax));
    
    
    
    if length(c_imgTrax) >0 %Check if valid tracks exist in this frame
        
        
        %% List XY Coords:
        for ctr2 = 1:length(c_imgTrax)
            tempIdx      = find(trackList{c_imgTrax(ctr2),1}(:,1)==ctr1);
            c_xy(ctr2,:) = [trackList{c_imgTrax(ctr2),1}(tempIdx,:)  c_imgTrax(ctr2)];
            c_xy(ctr2,2) = c_xy(ctr2,2)/scalingFactor;
            c_xy(ctr2,3) = c_xy(ctr2,3)/scalingFactor;
            
            %Trackmate sometimes sends along zeros as XY coords. Correcting for this:
            c_xy(c_xy(ctr2,2)==0,2)=1;
            c_xy(c_xy(ctr2,3)==0,3)=1;
        end
        
        
        %% Debugging: Overwrite XY points raw image:
        debugger =0;
        if debugger==1
            temp_xy      = [c_xy(:,2:3) c_xy(:,5)];  %Passing only XY and TrackID
            im_d1        = size(c_img,1);
            im_d2        = size(c_img,2);
            idxlist      = sub2ind([im_d1 im_d2], temp_xy(:,2), temp_xy(:,1));
            lbl          = c_img;
            lbl(idxlist) = 400;
        end
        
        keyboard
        %% Watershed Thresholding:
        temp_xy = [c_xy(:,2:3) c_xy(:,5)];  %Passing only XY and TrackID
        try
            [lbl, badtrax{ctr1}] = locwater3(c_img, temp_xy);
        catch
            keyboard
        end
        
        
        
        %% Data Extraction:
        try
            [tvec_med, tvec_mod, tvec_trax] = tmcalc(c_img, lbl);
        catch
            keyboard
        end
        
        
        
        %% Assigning into data structure:
        for ctr3 = 1:length(tvec_trax)
            resvec_seg(tvec_trax(ctr3)).median(ctr1) = tvec_med(ctr3);
            resvec_seg(tvec_trax(ctr3)).mode(ctr1)   = tvec_mod(ctr3);
        end
        
        
        clear tvec_med tvec_mod tvec_trax temp_xy c_xy c_imgTrax c_imgTrax1
        
        %% Debugging:
        ctr1
%         ff = figure('visible', 'off')
%         hold on
%         %subplot(1,2,1), imagesc(lbl), subplot(1,2,2), imagesc(c_img)
%         imagesc(lbl.*double(c_img))
%         print(ff,[c_dir 'segimg' num2str(ctr1)], '-dpng','-r200')
%         close gcf
        
        disp(['badpts: ' num2str(badtrax{ctr1})])
        
    end %If there are any valid tracks in the frame
end %Frame loop


save([c_dir 'tm_segmentervars.mat'],  'scalingFactor', 'invertmode' ,'trackList', 'badtrax', 'file' ,'resvec_seg', '-v7');





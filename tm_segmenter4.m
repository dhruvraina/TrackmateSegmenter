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
addpath('/Users/draina/Documents/MATLAB/Dhruv''s Stuff/MPI_Progs/ExchangeTools')


%% Inputs:
c_dir = '/Users/draina/Desktop/SiR_Test/';
file = [c_dir 'RegularOutput.xlsx'];

% c_dir = '/Users/draina/Desktop/Draq5_Test/';
% file = [c_dir 'InvertedTrack_matlab.xml'];
% [trackList trackMeta] = importTrackMateTracks(file);

rawTrackImport = readtable(file);

%Add +1 to the track_ID to deal with track_ID==0; This way I can uniformly
%subtract 1 from all track IDs.
rawTrackImport.TRACK_ID = rawTrackImport.TRACK_ID+1;


mainImg = tiffread2([c_dir 'SiR_Test_Nucleus.tif']);
measurementImg = tiffread2([c_dir 'SiR_Test_Sensor_inv.tif']);
scalingFactor = 1;  %:: Scaling factor is derived from image Metadata. TBD: automatically extract this number.


%% For Inverted Tracks:
% invertmode = 'y';
% framestart = 1;
% frameend = length(mainImg);
% if invertmode =='y'
%     
%     %Make a backup of the original variable for debugging
%     %templist = trackList;
%     
%     %Convert the frame numbers into reverse order
%     gg = cellfun(@(x) frameend-x(:,1), trackList, 'UniformOutput', 0);
%     
%     %Replace reversed frame numbers in original trackList
%     for cc = 1:length(trackList)
%         trackList{cc,1}(:,1)= gg{cc,1}(:,1);
%     end
%     
%     %Re-sort accordind to new frame numbers:
%     trackList = cellfun(@(x) sortrows(x,1), trackList, 'UniformOutput', 0);
%     
% end



%% Framewise Segmentation and Datextraction:
for ctr1 = 1:length(mainImg)
%for ctr1 = 140:230
  % ctr1 = 201;
    
    %Find all tracks which are active in this frame:
    c_img          = mainImg(ctr1).data;
    c_measureImg   = measurementImg(ctr1).data;
    c_imgTrax(:,1) = rawTrackImport.POSITION_X(rawTrackImport.FRAME==ctr1);
    c_imgTrax(:,2) = rawTrackImport.POSITION_Y(rawTrackImport.FRAME==ctr1);
    c_imgTrax(:,3) = rawTrackImport.TRACK_ID(rawTrackImport.FRAME==ctr1);
    
    if ~isempty(c_imgTrax) %Check if valid tracks exist in this frame
        
        %Check for weird zeros in XY Coords
        if any(c_imgTrax(:,1)==0) || any(c_imgTrax(:,2)==0)
            disp('Check XY Numbers for Zero Coordinates')
            keyboard
        end
        
        
        %% Debugging: Overwrite XY points raw image:
        debugger =0;
        if debugger==1
            im_d1        = size(c_img,1);
            im_d2        = size(c_img,2);
            idxlist      = sub2ind([im_d1 im_d2], c_imgTrax(:,1), c_imgTrax(:,2));
            lbl          = c_img;
            lbl(idxlist) = 400;
        end
        
        
        %% Nuclear Segmentation:
        try
            [lbl, badtrax{ctr1}] = nucseg(c_img, c_imgTrax);
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
        
        
        
        %% Debugging:
        ctr1
        ff = figure('visible', 'off');
        lbl2=lbl;
        lbl2(lbl==1000) = 0.2;
        lbl2(lbl2==0)=0.5;
        %lbl2(lbl2>-1) = 1;
        %subplot(1,2,1), imagesc(lbl), subplot(1,2,2), imagesc(c_img)
        %imshowpair(c_img, lbl2, 'ColorChannels', 'red-cyan')
        imagesc(lbl2)
        print(ff,[c_dir 'ImgSeg2/segimg2' num2str(ctr1)], '-dpng','-r200')
        close gcf
        
        disp(['badpts: ' num2str(badtrax{ctr1})]);

        clear tvec_med tvec_mod tvec_trax temp_xy c_xy c_imgTrax 

        
    end %If there are any valid tracks in the frame
    clear c_imgTrax
end %Frame loop


save([c_dir 'tm_segmentervars.mat'],  'scalingFactor', 'invertmode' ,'trackList', 'badtrax', 'file' ,'resvec_seg', '-v7');


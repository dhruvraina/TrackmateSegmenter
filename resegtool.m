function varargout = resegtool(im_d1, im_d2, xyimg, resegmenter, im, lblfinal, idxlist)

xyimg_dil = imdilate(xyimg, strel('disk', 5));
trackList = xyimg(xyimg>0);
%Resegmenter:
reseg_temp = zeros(im_d1, im_d2);

%Erase original mask, calculate Euler number
for cc = 1:length(resegmenter)
   
    breakup = resegmenter(cc);
    reseg_temp(lblfinal==breakup) = 1; %
    
   %Get trackIDs that need to be resegmented
    currxy = nonzeros(xyimg(reseg_temp==1));
       
    % Delete the mask for object to be resegmented (gets rid of all trackIDs in the masked region!)
    lblfinal(ismember(lblfinal,currxy)) = 0; 

    %Start the while loop
    rp_seg = reseg_temp;
    % masked_im = reseg_temp.*im;  %This doesn't need to be here. Its only there because I was thinking about implementing an 'intensity' based resegmenter in the while loop, instead of the current morphological one.
    [lbl num] = bwlabel(rp_seg);
    
    %Checking if trackIDs(in xyimg) are on different objects(lbl)
    xychecker = lbl(ismember(xyimg, currxy));
    dilator = 3;
    
    while length(unique(xychecker))~=length(currxy) %While trackIDs are on same obj.
        rp_seg = imerode(rp_seg, strel('disk', 1));
        rp_seg(ismember(xyimg_dil, currxy)) = 1;        %Forcing the actual XY pts to ALWAYS be there and avoid being eroded
        [lbl num] = bwlabel(rp_seg);                %bwlabel
        lbl(lbl>0) = lbl(lbl>0)+50;                 %Adding arbitrarily high number
        xychecker = lbl(ismember(xyimg, currxy));
        dilator = dilator+1;
        %figure, imagesc(lbl)
    end
    
    %bwareaopen returns a logical, I need to preserve lbl nums
    %lblclean = bwareaopen(lbl, 5);
    %lbl = lbl.*lblclean;

    
    lbldil = imdilate(lbl, strel('disk', dilator));
   
    %Reassign trackIDs to segmented objects
    for i = 1:length(xychecker)
    lbldil(lbldil==xychecker(i)) = nonzeros(xyimg(lbl==xychecker(i)));       %lbldil sometimes doesn't overlap with xypts due to mask alteration from imdilate. This is why I need the lbl image.
    end

    %Do the fucking addition. This is a shitty way of coding this because
    %it sometimes leads to carry-over numbers between loops.
    lblfin_t = zeros(im_d1, im_d2);
    lblfin_t(ismember(lbldil, xyimg)) = lbldil(ismember(lbldil, xyimg));
    lblfinal = lblfin_t+lblfinal;
    
    %Clear so this can run through the loop again
    reseg_temp = zeros(im_d1, im_d2);    
end

disp(['Resegmenting frame'])

%Delete all overlapping boundaries:
lblfinal(~ismember(lblfinal, trackList)) = 0;
varargout{1} = lblfinal;
end
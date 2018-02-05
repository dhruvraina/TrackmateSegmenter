%local watershed

function varargout = locwater(image, locxy)
generic = 0; %If generic on, it can do something else with outputs.

im = mat2gray(image);
locxy = round(locxy);
im = imtophat(im, strel('disk', 60));                                        %correct uneven background

im_d1 = size(im,1);
im_d2 = size(im,2);

%Get XY Points for tracks in image:
idxlist = sub2ind([im_d1 im_d2], locxy(:,2), locxy(:,1));
xypts = zeros(im_d1,im_d2);
xypts(idxlist) = 1;
xypts = imdilate(xypts, strel('disk', 2));


keyboard

%Create Masks based on XY Coords:
im_iobr = imreconstruct(xypts, im);
im_thresh = im2bw(im_iobr, graythresh(im_iobr)*1.5);
im_fil = imfill(im_thresh, 'holes');


%Running this twice in succession removes thinly connected objects
im_mask = imerode(im_fil, strel('disk', 2));
im_mask = imerode(im_mask, strel('disk', 1));
im_mask = imerode(im_mask, strel('disk', 1));
im_mask2 = imdilate(im_mask, strel('disk', 6));


%Watershedding:
im_mask3 = im_mask2.*im;
dist = bwdist(~im_mask3);                                                    %Distance transform for border lines
impmin = imimposemin(dist, im_mask | bwperim(im_mask2));                     %Impose min on borders and foreground
ww = watershed(impmin);
%figure, imagesc(ww)


%Catch XY coords just outside the segments:
ww_dil = imdilate(ww, strel('disk', 3));

%List the intersection points:
intersct = ww_dil(idxlist);


%Resegment if two xypts are in the same watershed basin:
resegmenter = sort(intersct);
resegmenter = double(resegmenter(diff(resegmenter)==0));


%If more than three spots in a mask: Error out
if ~isempty(double(resegmenter(diff(resegmenter)==0)))  %If more than two spots accidentally fall into the same bin
    keyboard
end

%If two spots in a mask:
if ~isempty(resegmenter)
    reseg_temp = zeros(im_d1, im_d2);
    
    %Erase original mask, calcular Euler number
    for cc = 1:length(resegmenter)
        breakup = resegmenter(cc);
        reseg_temp(ww==breakup) = 1; %
        ww_dil(ww_dil==breakup) = 0; %
        rp_seg = regionprops(reseg_temp, 'Euler');
        
        
        %Erode until change in Euler. Dilate later according to amount of erosion
        dilator = 3;   %dilator starts at ww_dil dilation strel size, and then grows
        while rp_seg.EulerNumber==1
            reseg_temp = imerode(reseg_temp, strel('disk', 1));
            rp_seg = regionprops(reseg_temp, 'Euler')
            dilator = dilator+1
        end
        
        %Clear small objects, relabel, dilate back to original size
        reseg_temp = bwareaopen(reseg_temp,5);
        [lbl num] = bwlabel(reseg_temp);
        lbl = imdilate(lbl, strel('disk',dilator));
        lbl = lbl*50;                 %arbitrarily high number prevents clash with existing indicies in Watershed img.
        ww_dil = double(ww_dil)+lbl;
        
        
        %Clear so this can run through the loop again
        reseg_temp = zeros(im_d1, im_d2);
    end
    
    %Update intersct
    intersct = ww_dil(idxlist);
    disp(['Resegmenting frame'])
    
end


%Area Filter
ww_rp = regionprops(ww_dil, 'Area');
area_filt = median([ww_rp(nonzeros(intersct)).Area]);                          %Area median of intersection between watershed label and xy points list
arfilt_H = area_filt*3;
arfilt_L = area_filt/4;

if generic ==0
    
    %Finding bad and good points:
    badpts_mask_all = find(([ww_rp(:).Area]>arfilt_H) | ([ww_rp(:).Area]<arfilt_L));
    badpts_locxy_idx = ismember(intersct, badpts_mask_all);
    %badpts3 = intersct(badpts2==1);                                             %Will I need this? Returns the idx of badpts in watershed img
    goodpts_locxy_all = intersct(badpts_locxy_idx==0);
    goodpts_locxy_idx = double(intersct).*~badpts_locxy_idx;
    locxy(:,4) = goodpts_locxy_idx;
        
    
    %Outputs:
    ww_dil(~ismember(ww_dil, goodpts_locxy_all)) = 0;
    goodtrax = nonzeros(~badpts_locxy_idx.*locxy(:,3));
    badtrax = nonzeros(badpts_locxy_idx.*locxy(:,3));

    %Changing lbls to reflect TrackIDs
    %cheating and using a for loop:
    ww_dil_cpy = ww_dil;
   % ww_dil = ww_dil_cpy;
    for aa = 1:size(locxy,1)
        ww_dil(ww_dil_cpy==locxy(aa,4)) = locxy(aa,3);
    end
    
end


varargout{1} = ww_dil;
varargout{2} = badtrax;

end
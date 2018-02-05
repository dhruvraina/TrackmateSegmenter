



function varargout = tmcalc(image, mask)

trackList = unique(nonzeros(mask(:)));

for aa = 1:length(trackList)
    resvec_med(aa,:) = median(image(mask==trackList(aa)));  %MEDIAN
    resvec_mod(aa,:) = mode(image(mask==trackList(aa)));    %MODE
    resvec_lis(aa,:) = trackList(aa);                       %TRACKLIST
end

varargout{1} = resvec_med;
varargout{2} = resvec_mod;
varargout{3} = resvec_lis;

end
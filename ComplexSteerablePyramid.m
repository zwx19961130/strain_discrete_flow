% function ComplexSteerablePyramid(I1, I2)
I1 = imread("data\test_warping1_discon1.jpg");
I2 = imread("data\test_warping2_discon1_predis16.5.jpg");

nF = 2;
refFrame = 1;

[h,w] = size(I1);


ht = maxSCFpyrHt(zeros(h,w));

filters = getFilters([h w], 2.^[0:-1:-ht], 4);

[croppedFilters, filtIDX] = getFilterIDX(filters);





reconLevel = @(im_dft, k) 2*(croppedFilters{k}.*fftshift(fft2(im_dft)));


numLevels = numel(filters);        
fprintf('Moving video to Fourier domain\n');
vidFFT = zeros(h,w,2);
% for k = 1:nF
%     originalFrame = rgb2ntsc(vid(:,:,k));
%     tVid = imresize(originalFrame(:,:,1), [h w]);
%     vidFFT(:,:,k) = single(fftshift(fft2(tVid)));
% end
% clear vid;

vidFFT(:, :, 1) = fftshift(fft2(I1));
vidFFT(:, :, 2) = fftshift(fft2(I2));

for level = 2:numLevels-1
    %% Compute phases of level
    % We assume that the video is mostly static
    pyrRef = buildLevel(vidFFT(:,:,refFrame), level);        
    pyrRef = angle(pyrRef);        

    delta = zeros(size(pyrRef,1), size(pyrRef,2) ,nF);
    fprintf('Processing level %d of %d\n', level, numLevels);
       
    
    filterResponse = buildLevel(vidFFT(:,:,2), level);
    pyrCurrent = angle(filterResponse);

    delta(:,:,frameIDX) = mod(pi+pyrCurrent-pyrRef,2*pi)-pi;                          
end


function result = buildLevel(im_dft, k)
    result = ifft2(ifftshift(croppedFilters{k}.* im_dft(filtIDX{k,1}, filtIDX{k,2})));
end


% end
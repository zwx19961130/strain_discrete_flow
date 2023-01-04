% function [u, v] = batch_demo(my_ksize, my_theta, my_gamma, my_sigma, my_gabor_lambda, my_psi, lambda, beta, my_pyramid_levels, my_pyramid_factor)
function RunTime = batch_demo_time(theta, w, ...
    lambda, beta, my_pyramid_levels, my_pyramid_factor)
    % load images
    %fpath = './frame10.png';
    %fpath = './Mars-1.jpg';
    addpath(genpath('.'));
    
    % fpath = './data/0275.png';
    fpath = "./data/test_warping1_discon1.jpg";
    img_src1 = im2double(imread(fpath));
    if size(img_src1,3) == 3
        img_src1 = rgb2gray(img_src1);
    end
    %fath = './frame11.png';
    %fpath = './Mars-2.jpg';
    % fpath = './data/0279.png';
    fpath = "./data/test_warping2_discon1_predis16.5.jpg";
    img_src2 = im2double(imread(fpath));
    if size(img_src2,3) == 3
        img_src2 = rgb2gray(img_src2);
    end
    %figure; imshow(img_src1);
    %figure; imshow(img_src2);
    
    % tv-l1 flow (coarse-to-fine)
    options.lambda = lambda;
    options.beta   = beta;
    max_iter = 50;
    options.max_iter = round(max_iter);
    check = 10;
    options.check = round(check);
    pyramid_levels = my_pyramid_levels;
    options.pyramid_levels = round(pyramid_levels);
    options.pyramid_factor = my_pyramid_factor;
    warps = 1;
    options.warps = round(warps);
    
    
%     options.gabor.ksize = my_ksize;
%     options.gabor.theta = my_theta;
%     options.gabor.gamma = my_gamma;
%     options.gabor.sigma = my_sigma;
%     options.gabor.gabor_lambda = my_gabor_lambda;
%     options.gabor.psi = my_psi;

% 
    options.steerable.theta = theta;
    options.steerable.w = w;

     


% x1 = [0.00940000000000000	0.114800000000000	0.396400000000000	0.0601000000000000	0.921300000000000	0.0601000000000000	0.396400000000000	0.114800000000000	0.00940000000000000]';
% x2 = [0.000800000000000000	0.0176000000000000	0.166000000000000	0.638300000000000	1	0.638300000000000	0.166000000000000	0.0480000000000000	0.000800000000000000]';
% y2 = [0.000800000000000000	0.0176000000000000	0.166000000000000	0.638300000000000	1	0.638300000000000	0.166000000000000	0.0176000000000000	0.000800000000000000]';
% y1 = [-0.00980000000000000	-0.0618000000000000	0.0998000000000000	0.755100000000000	0	-0.755100000000000	-0.0998000000000000	0.0618000000000000	0.00980000000000000]';
% R_h = x2 * x1.';
% I_h = y2 * y1.';
% 
% 
% I1 = img_src1;
% I2 = img_src2;
% 
% image1_horizontal = conv2(I1, R_h, 'same') + 1i * conv2(I1, I_h, 'same');
% image1_phase_horizontal = angle(image1_horizontal);
% 
% 
% 
% image2_horizontal = conv2(I2, R_h, 'same') + 1i * conv2(I2, I_h, 'same');
% image2_phase_horizontal = angle(image2_horizontal);
% 
% img_src1 = image1_phase_horizontal;
% img_src2 = image2_phase_horizontal;
    
    f = @() batch_tv_l1_optical_flow_time(img_src1, img_src2, options);
    
    RunTime = timeit(f);
    
    
    % post-processing
    % find robust max flow for better visualization
%     magnitude = (u.^2 + v.^2).^0.5;  
%     max_flow = prctile(magnitude(:),95);
%     
%     illumination = illumination - min(illumination(:));
%     illumination = illumination/max(illumination(:));
%     
%     u = min(max(u,-max_flow),max_flow);
%     v = min(max(v,-max_flow),max_flow);
    
    %% write flow
    %motfile = 'motion.png';
    %illfile = 'illumination.png';
    %[M N C] = size(img_src1);
    %flow = zeros(M,N,2);
    %flow(:,:,1) = u;
    %flow(:,:,2) = v;
    %imwrite(uint8(flowToColor(flow)),motfile);
    %imwrite(illumination,illfile);
    
    % append the estimated optical flow to the first frame to show the motion
    % of every pixel


%     
%     figure; imshow(img_src1);
%     hold on;
%     [n_row,n_col] = size(img_src1);
%     [X,Y] = meshgrid(1:10:n_col,1:10:n_row);
%     quiver(X,Y,u(1:10:end,1:10:end),v(1:10:end,1:10:end));
%     title('estimated optical flow appended to the first frame');
    
    % warping
    %sz = size(img_src1);
    %[X,Y] = meshgrid(1:sz(2), 1:sz(1));
    %back = interp2(img_src2,X+flow(:,:,1),Y+flow(:,:,2));
    %figure; imagesc(abs(back-img_src1))



%     im1_back = warpImage(img_src2, u, v);
%     figure; imshow(im1_back); title('warped I2 -> I1');
end

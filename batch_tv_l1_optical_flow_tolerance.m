function [u, v, w, p] = batch_tv_l1_optical_flow_tolerance(I1, I2, options) 
% Inputs: I1, I2 are two consecutive images (only grayscale allowed)
%         options.lambda (10-100)
%         options.beta (0.001,0.01,0.1)
%         options.warps (1-5)
%         options.max_iter (50-200)
%         options.pyramid_levels (100-1000)
%         options.pyramid_factors (0.5-0.9)
% Outputs: (u,v) is the flow, where size(u) = size(v) = size(I1)
%          w is the estimated illumination changes

lambda = options.lambda;
beta = options.beta;
warps = options.warps;
%max_iter = options.max_iter;
pyramid_levels = options.pyramid_levels;
pyramid_factor = options.pyramid_factor;
check = options.check;
tolerance = options.tolerance;

%k=1;

[M,N,C] = size(I1);

num_dual_vars = 6;

width_Pyramid = cell(pyramid_levels,1);
height_Pyramid = cell(pyramid_levels,1);

I1_Pyramid = cell(pyramid_levels,C);
I2_Pyramid = cell(pyramid_levels,C);

% precompute image sizes
width_Pyramid{1} = N;
height_Pyramid{1} = M;
for i = 2:pyramid_levels
  width_Pyramid{i} = pyramid_factor*width_Pyramid{i-1};
  height_Pyramid{i} = pyramid_factor*height_Pyramid{i-1};
  tolerance = tolerance*2;
  if min(width_Pyramid{i}, height_Pyramid{i}) < 16
    pyramid_levels = i;
    break;
   end
end
for i = 1:pyramid_levels
  width_Pyramid{i} = round(width_Pyramid{i});
  height_Pyramid{i} = round(height_Pyramid{i});
end

% set up image pyramides
for i = 1:pyramid_levels
  if i == 1
    for j=1:C
      I1_Pyramid{1,j} = I1(:,:,j);
      I2_Pyramid{1,j} = I2(:,:,j);

      I1_Pyramid{1,j} = fftshift(fft2(I1_Pyramid{1,j}));
      I2_Pyramid{1,j} = fftshift(fft2(I2_Pyramid{1,j}));
    end
  else
    for j=1:C
      I1_Pyramid{i,j} = imresize(I1_Pyramid{i-1,j}, ...
                                 [height_Pyramid{i} width_Pyramid{i}], 'bicubic');
      I2_Pyramid{i,j} = imresize(I2_Pyramid{i-1,j}, ...
                                 [height_Pyramid{i} width_Pyramid{i}], 'bicubic'); 

      I1_Pyramid{1,j} = fftshift(fft2(I1_Pyramid{1,j}));
      I2_Pyramid{1,j} = fftshift(fft2(I2_Pyramid{1,j}));
    end
  end
end

% main loop
for level = pyramid_levels:-1:1  %1--->2
  %scale = pyramid_factor^(level-1);
  
  M = height_Pyramid{level};
  N = width_Pyramid{level};
 
  if level == pyramid_levels
    % initialization  
    u = zeros(M,N);
    v = zeros(M,N);
    w = zeros(M,N);
    p = zeros(M,N,num_dual_vars);   
  else
    rescale_factor_u = width_Pyramid{level+1}/width_Pyramid{level};
    rescale_factor_v = height_Pyramid{level+1}/height_Pyramid{level};
    
    % prolongate to finer grid    
    u = imresize(u,[M N], 'nearest')/rescale_factor_u;    
    v = imresize(v,[M N], 'nearest')/rescale_factor_v;
    w = imresize(w,[M N], 'nearest');
    
    p_tmp = p;
    p = zeros(M,N,num_dual_vars); 
    for i=1:num_dual_vars
      p(:,:,i) = imresize(p_tmp(:,:,i),[M N], 'nearest');
    end
  end
  
  I1 = zeros(M,N,C);
  I2 = zeros(M,N,C);
  for j=1:C
    I1(:,:,j) = I1_Pyramid{level,j};
    I2(:,:,j) = I2_Pyramid{level,j};
  end
  
  fprintf('*** level = %d\n', level);



[J1,~,~,~,~,~,~] = steerGaussianFilter(I1, options.steerable.theta, options.steerable.w, 0);
[J2,~,~,~,~,~,~] = steerGaussianFilter(I2, options.steerable.theta, options.steerable.w, 0);


image1_phase_horizontal = angle(J1);
image2_phase_horizontal = angle(J2);


  
  [u, v, w, p, k] = champolle_pock_primal_dual_alg(image1_phase_horizontal, image2_phase_horizontal, ...
                        u, v, w, p, lambda, warps, beta, tolerance, check); 
                    %check, max_iter, scale);
 
  fprintf('The number of iteration steps is %d\n', k);
                    
  tolerance = tolerance/2;
  
  fprintf('=== Processing...%.2f/100\n', (pyramid_levels + 1 - level)/pyramid_levels*100);
  drawnow;
end
fprintf('Done!\n');
drawnow;
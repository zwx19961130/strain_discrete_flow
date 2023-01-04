
% gabor_lambda = 0.3 : 0.3 : 3;

% ksize = 2;
% theta = 121;
% gamma = -12;
% sigma = 1.8;
% gabor_lambda = 1.5;
% psi = 4.9;


theta = 23;
w = 1.23;




% w = 3 : 0.1 : 6;
lambda = 20;
beta = 0.01;
pyramid_level = 6;
pyramid_factor = [0.2, 0.5, 1];


for i1 = 1 : length(pyramid_factor)
%         for i2 = 1 : length(w)
%            for i4 = 1 : length(beta)
%                for i5 = 1 : length(pyramid_level)
%                    for i6 = 1 : length(pyramid_factor)


%     [u, v] = batch_demo(ksize, deg2rad(theta),...
%         gamma, sigma, gabor_lambda, psi,...
%         lambda, beta, pyramid_level, pyramid_factor);


RunTime = batch_demo_time(theta, w, ...
    lambda, beta, pyramid_level, pyramid_factor(i1));

% close all

%  sum = 0;
%  for i = 1 : 500
%      for j = 1 : 500
%         sum = sum + abs(sqrt(u_ROI(i, j) ^ 2 + v_ROI(i, j) ^ 2) - sqrt(u_hat(i, j) ^ 2 + v_hat(i, j) ^ 2));
%      end
%  end



%  MDE = sum(abs(sqrt(u_ROI .^ 2 + v_ROI .^ 2) - sqrt(u_hat .^ 2 + v_hat .^ 2)), "all") / (500 * 500);
 
 
 fileID =fopen("batch_test.txt", "a");
%  
%      formatSpec = "%d\t%d\t%f\t%f\t%f\t%f\t" + ...
%          "%d\t%f\t%d\t%f\t" + ...
%          "%f\t%f\n";
% 
%  fprintf(fileID, formatSpec, ksize, theta(i1), gamma, sigma, gabor_lambda, psi, lambda, beta, pyramid_level, pyramid_factor, rmse, MDE);

formatSpec = "%d\t%f\t" + ...
    "%d\t%f\t%d\t%f\t%f\n";
 fprintf(fileID, formatSpec, theta, w, lambda, beta, pyramid_level, pyramid_factor(i1), RunTime);
 
%                    end
%                end
%            end
%         end
fclose(fileID);
%             end
%         end
%     end
end








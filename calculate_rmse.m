
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
pyramid_factor = 0.5;


for i1 = 1 : length(lambda)
%         for i2 = 1 : length(w)
%            for i4 = 1 : length(beta)
%                for i5 = 1 : length(pyramid_level)
%                    for i6 = 1 : length(pyramid_factor)


%     [u, v] = batch_demo(ksize, deg2rad(theta),...
%         gamma, sigma, gabor_lambda, psi,...
%         lambda, beta, pyramid_level, pyramid_factor);

[u, v] = batch_demo(theta, w, ...
    lambda(i1), beta, pyramid_level, pyramid_factor);
u_ROI = u(21:520, 21:520); % 取出Region of Interest
v_ROI = v(21:520, 21:520);


u_hat = zeros(500, 500); % 新建u_hat矩阵，u_hat表示ground truth值
u_hat(:, 250:end) =0.5;
v_hat = u_hat;


numerator = sum((u_hat - u_ROI) .^ 2, 'all') + sum((v_hat - v_ROI) .^ 2, 'all'); % 分子部分，把u和v的2次方误差加起来
rmse = sqrt(numerator / (500 * 500)); % u的大小是500 * 500, 因为u,v 2个矩阵，所以乘以2
close all

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
 fprintf(fileID, formatSpec, theta, w, lambda(i1), beta, pyramid_level, pyramid_factor, rmse);
 
%                    end
%                end
%            end
%         end
fclose(fileID);
%             end
%         end
%     end
end








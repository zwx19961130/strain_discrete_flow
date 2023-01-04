function J = steerGaussFilterOrder2_Boualem(I,theta,w,disFlag,sigma,choice)

%    J = steerGaussFilterOrder2(I,theta,sigma,disFlag) 
%    calculated 
%    Input:
%    1. I: input image
%    2. theta: the orientation
%    3. sigma: standard deviation of the Gaussian template   
%    Output:
%    J. The response of derivative in theta direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if choice ==1  %%% 1st Version

if(mod(w,2)~=0)
        'error';
end
    x = [-w:w];
    [xx,yy] = meshgrid(x,x);
    sx = 1/sqrt(0.9213); sy = sx;
    G2_0 = ( ((1/4)/pi) /sx^6 )*(2*xx.^2-2*sx^2).*exp(-0.5*((xx./sx).^2 + (yy/sy).^2));   % Freman paper
    G2_90 =    1/4/pi   /sy^6  *(2*yy.^2-2*sy^2).*exp(-0.5*((xx./sx).^2 + (yy/sy).^2));
%%%% Normalization
    %fun = @(x0,y0) (1/4/pi/sx^6*(2*x0.^2-2*sx^2).*exp(-0.5*((x0./sx).^2 + (y0/sy).^2))).^2;
    %q = integral2(fun,-inf,inf,-inf,inf);
    %s = sqrt(1/q);
    %normalisation
    %G2 = G2*s;
%%%%   Hilbert  
    for i = 1:length(x)
        H2_0(i,:) = imag(hilbert(G2_0(i,:)));
        H2_90(:,i) = imag(hilbert(G2_90(:,i)));
    end
%%%%   Filtering
    Ix_g = imfilter(I,G2_0,'same','replicate');
    Iy_g = imfilter(I,G2_90,'same','replicate');
    Ix_h = imfilter(I,H2_0,'same','replicate');
    Iy_h = imfilter(I,H2_90,'same','replicate');
    
    Jg = cos(theta)*Ix_g+sin(theta)*Iy_g;
    Jh = cos(theta)*Ix_h+sin(theta)*Iy_h;
    J = Jg + 1i*Jh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else %% 2st Version
%w=floor((8/2)*sigma);     
% if w < 1
%   w = 1;
% end
x = [-3.3500 -2.6800 -2.0100   -1.3400   -0.6700        0    0.6700    1.3400    2.0100  2.6800 3.3500];%2.0100
[xx,yy] = meshgrid(x,x);

 G2b = exp(-(xx.^2+yy.^2)).*(0.9213*(2.*yy.^2-1));
 G2a = (exp(-(xx.^2+yy.^2)).*(0.9213*(2.*xx.^2-1)))/0.9780;
 H2b =exp(-(xx.^2+yy.^2)).*(-2.254.*yy+yy.^3).*0.9780;    %exp(-(xx.^2+yy.^2)).*(-.7515+xx.^2).*yy.*0.9780;%
 H2a =exp(-(xx.^2+yy.^2)).*(-2.254.*xx+xx.^3).*0.9780;    %G2a;%;

%      for i = 1:length(x)
%         H2a(i,:) = imag(hilbert(G2a(i,:)));
%         H2b(:,i) = imag(hilbert(G2b(:,i)));
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: Determine oriented filter response.
% Calculate image gradients (using separability).
Ix_g = imfilter(I,G2a,'same','replicate');
Iy_g = imfilter(I,G2b,'same','replicate');
Ix_h = imfilter(I,H2a,'same','replicate');
Iy_h = imfilter(I,H2b,'same','replicate');
% Evaluate oriented filter response.
%J = (cos(theta))^2*I2a+sin(theta)^2*I2c-2*cos(theta)*sin(theta)*I2b;
    Jg = cos(theta)*Ix_g+sin(theta)*Iy_g;
    Jh = cos(theta)*Ix_h+sin(theta)*Iy_h;
    J = Jg + 1i*Jh;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%     if disFlag
%         figure(1); 
%         clf; 
%         set(gcf,'Name','Oriented Filtering');
%         subplot(2,3,1); imagesc(I);  colormap(gray);
%         title('Input Image');
%         subplot(2,3,2); imagesc(abs(J)); colormap(gray);
% %         hold on
% %         scatter(x0(1:15),y0(1:15),'r+','linewidth',3);
%         title(['Filtered Image Amplitude(\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
%         subplot(2,3,3); 
%         imagesc(angle(J));axis image; colormap(gray);
%         title(['Filtered Image Phase (\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
%         subplot(2,3,4); 
%         imagesc(G2_0);axis image; colormap(gray);
%         title(['Oriented Filter G2 (\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
%         subplot(2,3,5); 
%         imagesc(H2_0);axis image; colormap(gray);
%         title(['Oriented Filter H2 (\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
%     end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
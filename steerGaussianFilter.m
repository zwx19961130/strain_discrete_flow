% Implementation of steerable filter 
function [J,Jg,Jh,G2_0,G2_90,H2_0,H2_90] = steerGaussianFilter(I,theta,w,disFlag)
if(mod(w,2)~=0)
        'error';
end
    x = [-w:w]*0.67;
    [xx,yy] = meshgrid(x,x);
    % Freeman's equation
    G2_0 = 0.9213*(2*xx.^2-1).*exp(-xx.^2-yy.^2);
    G2_90 = 0.9213*(2*yy.^2-1).*exp(-xx.^2-yy.^2);
    H2_0 = (-2.205*xx+0.978*xx.^3).*exp(-xx.^2-yy.^2);
    H2_90 = (-2.205*yy+0.978*yy.^3).*exp(-xx.^2-yy.^2);
    
%    for i = 1:length(x)
%         H2_0(i,:) = imag(hilbert(G2_0(i,:)));
%         H2_90(:,i) = imag(hilbert(G2_90(:,i)));
%     end
    Ix_g = imfilter(I,G2_0,'same','replicate');
    Iy_g = imfilter(I,G2_90,'same','replicate');
    Ix_h = imfilter(I,H2_0,'same','replicate');
    Iy_h = imfilter(I,H2_90,'same','replicate');
    
    
    
    
    Jg = cos(theta)*Ix_g+sin(theta)*Iy_g;
    Jh = cos(theta)*Ix_h+sin(theta)*Iy_h;
    J = Jg + 1i * Jh;
    
    if disFlag
        figure(1); 
        clf; 
        set(gcf,'Name','Oriented Filtering');
        subplot(2,3,1); imagesc(I);  colormap(gray);
        title('Input Image');
        subplot(2,3,2); imagesc(abs(J)); colormap(gray);
%         hold on
%         scatter(x0(1:15),y0(1:15),'r+','linewidth',3);
        title(['Filtered Image Amplitude(\theta = ',num2str(theta*(180/pi)),'{\circ})']);
        subplot(2,3,3); 
        imagesc(angle(J)); colormap(gray);
        title(['Filtered Image Phase (\theta = ',num2str(theta*(180/pi)),'{\circ})']);
        subplot(2,3,4); 
        imagesc(G2_0); colormap(gray);
        title(['Oriented Filter G2 (\theta = ',num2str(theta*(180/pi)),'{\circ})']);
        subplot(2,3,5); 
        imagesc(H2_0); colormap(gray);
        title(['Oriented Filter H2 (\theta = ',num2str(theta*(180/pi)),'{\circ})']);
    end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
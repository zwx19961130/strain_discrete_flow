function gb = my_gabor(ksize, theta, gamma, sigma, lambda, psi)
    d = (ksize + 1) / 2;
    gb = zeros(ksize, ksize);
    for x = 1 : ksize
        xd = x - d; % distance from the center
        for y = 1 : ksize
            yd = y - d;
            xn = xd*cos(theta) + yd*sin(theta);
            yn = -xd*sin(theta)+ yd*cos(theta);
            gb(x, y)= exp(-(xn^2+gamma^2*yn^2)/(2*sigma^2))*exp(1i*(2*pi*xn/lambda + psi));
        end
    end
end
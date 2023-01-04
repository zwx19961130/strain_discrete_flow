theta = 0;
gamma = 1;
psi = 1;
sigma = 1;
lambda = 1;




gb = gabor(theta, gamma, sigma, lambda, psi);

imagesc(real(gb))
colorbar
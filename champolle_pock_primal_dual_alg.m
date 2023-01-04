function [u, v, w, p, k] = champolle_pock_primal_dual_alg(I1, I2, u, v, w, p, lambda, warps, beta, tolerance, check)
                           %maxits, scale, check)
% [M N C] = size(I1);
[n_row,n_col,c] = size(I1);
N             = n_row * n_col;
nabla         = make_nabla(n_col,n_row); % Kxf, Kyf for gradient
num_dual_vars = 3;

% stepwidth
L = sqrt(8);  %Lipschitz constant of gradient

tau = 1/L;
sigma = 1/L;

% some parameters
epsilon_u = 0;
epsilon_w = 0;
% vectorization
vector_u   = reshape(u',N,1); % primal
vector_v   = reshape(v',N,1); % primal
vector_w   = reshape(w',N,1); % auxilary primal illumination
vector_u_  = vector_u; % primal_bar
vector_v_  = vector_v; % primal_bar
vector_w_  = vector_w; % primal_bar
vector_p   = zeros(2*N,num_dual_vars);
for i = 1:num_dual_vars
  vector_p(1:N,i)     = reshape(p(:,:,2*i-1)',N,1);
  vector_p(N+1:2*N,i) = reshape(p(:,:,2*i)',N,1);
end
nabla_t = nabla'; % K^T for divergence

% inner loops for warping
for j = 1:warps
  %vector_u0 = zeros(N,1);
  %vector_v0 = zeros(N,1);
  vector_u0 = vector_u;
  vector_v0 = vector_v;
  u0 = (reshape(vector_u0,n_col,n_row))';
  v0 = (reshape(vector_v0,n_col,n_row))';
  
  % warping
  fprintf('tv-l1-of-pd: warp = %d\n', j);
  
  [I_x, I_y, I_t, I2_warped] = warping(I1, I2, vector_u0, vector_v0);


  % vectorization
  vector_I_x = reshape(I_x',N,1);
  vector_I_y = reshape(I_y',N,1);
  vector_I_t = reshape(I_t',N,1);
  I_grad_sqr = max(1e-09, vector_I_x.^2 + vector_I_y.^2 + beta*beta);
  
 % iterations of Champolle-Pock alg
 %for k = 0:1:10
  k=1;
  while 1
    %% DUAL
    % compute derivatives and update dual variable: p_tilde = p + sigma*K*uv_bar
    vector_p(:,1:2) = (vector_p(:,1:2) + sigma*nabla*[vector_u_,vector_v_])/(1 + sigma*epsilon_u);
    vector_p(:,3)   = (vector_p(:,3)   + sigma*nabla*vector_w_)/(1 + sigma*epsilon_w); 
    
    % reprojection to |pu| <= 1
    norm = repmat(sqrt(vector_p(1:N,1).^2 + vector_p(N+1:2*N,1).^2 + ...
        vector_p(1:N,2).^2 + vector_p(N+1:2*N,2).^2),[2,2]);
    reprojection = max(1.0,norm);
    vector_p(:,1:2) = vector_p(:,1:2)./reprojection;
    
    norm = repmat(sqrt(vector_p(1:N,3).^2 + vector_p(N+1:2*N,3).^2),[2,1]);
    reprojection = max(1.0,norm);
    vector_p(:,3) = vector_p(:,3)./reprojection;
    
    %% PRIMAL
    % remember old u,v,w
    vector_u_ = vector_u;
    vector_v_ = vector_v;
    vector_w_ = vector_w;
  
    % uv_tilde = uv - tau*K^T*p
    vector_u = vector_u - tau * nabla_t * vector_p(:,1);
    vector_v = vector_v - tau * nabla_t * vector_p(:,2);
    vector_w = vector_w - tau * nabla_t * vector_p(:,3);
    
    % prox operator for u,v,w
    rho = vector_I_t + (vector_u - vector_u0).*vector_I_x + ...
        (vector_v - vector_v0).*vector_I_y + beta*vector_w;
    %I think u0,v0 and w0 is constant
    idx1 = rho      < - tau*lambda*I_grad_sqr;
    idx2 = rho      >   tau*lambda*I_grad_sqr;
    idx3 = abs(rho) <=  tau*lambda*I_grad_sqr;
    
    vector_u(idx1) = vector_u(idx1) + tau*lambda*vector_I_x(idx1);
    vector_v(idx1) = vector_v(idx1) + tau*lambda*vector_I_y(idx1);
    vector_w(idx1) = vector_w(idx1) + tau*lambda*beta;
    
    vector_u(idx2) = vector_u(idx2) - tau*lambda*vector_I_x(idx2);
    vector_v(idx2) = vector_v(idx2) - tau*lambda*vector_I_y(idx2);
    vector_w(idx2) = vector_w(idx2) - tau*lambda*beta;
    
    vector_u(idx3) = vector_u(idx3) - rho(idx3).*vector_I_x(idx3)./I_grad_sqr(idx3);
    vector_v(idx3) = vector_v(idx3) - rho(idx3).*vector_I_y(idx3)./I_grad_sqr(idx3);
    vector_w(idx3) = vector_w(idx3) - rho(idx3).*beta./I_grad_sqr(idx3);
    
    vector_u_ = 2*vector_u - vector_u_;
    vector_v_ = 2*vector_v - vector_v_;
    vector_w_ = 2*vector_w - vector_w_;
    
    
    if mod(k,check) == 0
     u = (reshape(vector_u,n_col,n_row))';
     v = (reshape(vector_v,n_col,n_row))';
     w = (reshape(vector_w,n_col,n_row))';      
     figure(109); show_flow(u,v,beta*w,I1,I2_warped + (u-u0).*I_x + (v-v0).*I_y + beta*w);
      fprintf('tv-l1-motion-primal-dual: it = %d\n', k)
    end
    
    tmp=zeros(N,3);
    tmp(:,1) = vector_u_-vector_u;
    tmp(:,2) = vector_v_-vector_v;
    tmp(:,3) = vector_w_-vector_w;
    
    a = tmp(:,1);
    b = tmp(:,2);
    c = tmp(:,3);
    
    s = a.^2+b.^2+c.^2;
    as = mean(s(:));
    
    if as <= tolerance^2
        break;
    end
    
    k=k+1;

  end
  
  % recover matrix representation
  u = (reshape(vector_u,n_col,n_row))';
  v = (reshape(vector_v,n_col,n_row))';
  w = (reshape(vector_w,n_col,n_row))';
  for i = 1:num_dual_vars
      p(:,:,2*i-1) = (reshape(vector_p(1:N,i),n_col,n_row))';
      p(:,:,2*i)   = (reshape(vector_p(N+1:2*N,i),n_col,n_row))';
  end
  % filter strong outliers
  u = peakfilt(u);
  v = peakfilt(v);
end

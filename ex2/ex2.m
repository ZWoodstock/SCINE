close all; clear all;
page_output_immediately(1)

load -mat 'xtrue_ex2.dat'
res = length(x_true);
nxt = norm(x_true);

%Proxified degradation:
gam_2 = .05;
F1 = @(x) prox_abs(x,gam_2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F2 subdivides x and approximates it with a step function
%based on subinterval averaging. num_subdivisions should be a power
%of 2, since it must divide the length of x_true (1024).
subdiv_num = 16;
subdiv_length = res/subdiv_num;
if length(x_true) / subdiv_num ~= subdiv_length
    fprintf('Improper selection of step function approximation!\n')
end
F2 = @(x) step_appr(x,subdiv_num,subdiv_length);
%Noise
noise2_param = 0.3;
noise2 = noise2_param*randn(size(x_true));
if noise2_param~=0
    snr2 = 20*log10(norm(x_true(:))./norm(noise2(:)))
end
p2 = F2(x_true+noise2);
T2 = @(x) x - F2(x) + p2;

gam_gradient = 0.025;
F3 = @(x) boundgrad(x,gam_gradient);
%p3 = 0;
T3 = @(x) x - F3(x);

%Generate data based on random inner products.
num_rand_obs = 1200;
%E holds the random vectors used to generate inner products with x_true
E = randn(length(x_true),num_rand_obs);
%rho holds the nonlinearly degraded inner products with the columns
%of E. Since one knows soft-thresholded rho if and only if they
%know the Tao-Vidakovic thresholded rho, we just store the
%proxified versions:

rho = zeros(num_rand_obs,1);
for i = 1:num_rand_obs
    E(:,i) = E(:,i)/norm(E(:,i));
    %Inner product step for rho:
    rho(i) = dot(E(:,i),x_true);
end
noise1_param = 0.025;
noise1 = noise1_param*randn(size(rho));
%Modeled Nonlinearity
NL = @(x) (abs(x)>gam_2).*(sign(x).*sqrt(x.^2-gam_2.^2));
%Proxifier for the modeled nonlinearity
S = @(x) sign(x).*(sqrt(x.^2+gam_2.^2)-gam_2);
%Nonlinearity which actually creates our data:
NL_real = @(x) (abs(x)>gam_2).*(sign(x).*nthroot(abs(x.^4-gam_2.^4),4));

rho = S(NL_real(rho)) + noise1;
if noise1_param ~= 0
    snr1 = 20*log10(norm(rho(:))./norm(noise1(:)))
end
%block_size is the number of firmly nonexpansive maps activated at
%each iteration (and is taken to be constant for this experiment).
block_size = num_rand_obs; %300;
num_blocks = num_rand_obs / block_size;
if num_blocks~=round(num_blocks)
    fprintf('ERROR: Inconsistent block structure')
end
%Initialize variables for loop
current_block = 0;
min_err = Inf;
min_x  = zeros(size(x_true));
x0 = zeros(size(x_true));
xold = x0;
eold = Inf;
err = Inf;
tol = 10^-3;
max_it = 50001/4;
k=0;
xold = zeros(size(x_true));
nn = 1/(2+num_rand_obs);

fprintf('\n Starting main. alg...\n')  

while k<max_it
    block_ind = block_size*current_block;
    for i = 1:block_size
        r_blocks_x(:,block_ind+i) = xold + E(:,block_ind+i).*...
           (rho(block_ind+i) - F1(dot(xold,E(:,block_ind+i))));
    end
    t2x = T2(xold);
    t3x = T3(xold);

    xnew = nn*(t2x+t3x+sum(r_blocks_x,2));
    xold=xnew;
    k=k+1;
    current_block = mod(current_block+1,num_blocks);
    if mod(k,500)==0
      fprintf('ITER: %2.1f %%\t ERR: %2.6f %%\n',100*k/max_it,...
                100*norm(xnew-x_true)/nxt)
    end
end

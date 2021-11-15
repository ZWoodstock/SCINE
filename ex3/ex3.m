close all; clear all;
pkg load image
page_output_immediately(1)

x_true = load('words256.mat');
nxt = norm(x_true,'fro');

sparsity_param = 1.5;
T1 = @(x) prox_abs(x,sparsity_param);
T1_skip = [1];

%Compression as in Andrews, '76. + blurring precomposition
tsvd_thr = 500;
windowsize2 = 7; %Windowsize should be an ODD number
ker2 = fspecial('average',windowsize2);
F2 = @(x) imfilter(tsvd(imfilter(x,ker2,'conv','replicate')...
           ,tsvd_thr,'soft'),ker2,'conv','replicate');
R2 = @(x) tsvd(imfilter(x,ker2,'conv','replicate'),tsvd_thr);
noise2_param = 12;
noise2 = noise2_param*randn(size(x_true));
x2_corrupted = sat(x_true+noise2,255,0);
noise2 = x_true - x2_corrupted;
if noise2_param~=0
    SNR2 = 20*log10(norm(x_true,'fro')./norm(noise2,'fro'))
end
[u, s, v] = svd(x_true+noise2);
tsvd_rank = numel(find(abs(diag(s)) >= tsvd_thr))
%T2 will be reactivated every 'T2_skip' iteration(s).
T2_skip = [5];
p2 = F2(x_true+noise2);
p2_show = R2(x_true+noise2);

%Initialize variables for loop

xold = zeros(size(x_true));
xold2=xold;
xold3=xold;
err = Inf;
tol = 10^-6;
k=0;
[N M]=size(xold);
R1x = zeros(N,M,2);
max_it = 1500;
E = zeros(max_it,1);
k=0;
T2 = @(x) p2 + x - F2(x);
T3 = @(x) sat(x,255,0);

fprintf('\n Starting main. alg...\n')  
R1x(:,:,1) = T1(xold);
R1x(:,:,2) = T2(xold);

nn2 = 1/2;
ds_coeff = 0.5;
coeff2 = 1-ds_coeff;
while k<max_it
    if mod(k,T1_skip(1))==0
      R1x(:,:,1) = T1(xold);
    end
    if mod(k,T2_skip(1))==0
      R1x(:,:,2) = T2(xold);
    end
    xnew = T3(ds_coeff*R1x(:,:,1) + (coeff2)*R1x(:,:,2));
    E(k+1) = norm(xnew-x_true,'fro');
    xold=xnew;
    k=k+1;
    if mod(k,150)==0
      fprintf('%2.1f %%\n',100*k/max_it)
    end
end
  
figure(1)
subplot(1,3,1)
imshow(uint8(x_true))
title('Original')

subplot(1,3,2)
imshow(uint8(p2_show))
title(strcat('Rank-',num2str(tsvd_rank),' Approximation'))

subplot(1,3,3)
imshow(uint8(xnew))
title('Recovery')

final_errs=100.*norm(xnew-x_true,'fro')./nxt;
fprintf(cstrcat('\n',num2str(max_it),' Iterations,\t', ...
      num2str(floor(max_it/T2_skip(1))),' SVDs,',...
      num2str(final_errs),'%% error\n'))







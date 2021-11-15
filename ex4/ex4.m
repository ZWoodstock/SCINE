clear all; close all;
pkg load image
pkg load signal
page_output_immediately(1)

X = rgb2gray(imread('NASA.tif'));
N = 600;
X = double(imresize(X,[N,N]));
N = length(X);
nxt = norm(X(:));

tsvd_thr = 1500;
F2 = @(x) tsvd(x(:,:,1)+x(:,:,2),tsvd_thr,'soft');
R2 = @(x) tsvd(x(:,:,1)+x(:,:,2),tsvd_thr);
noise2_param = 0;
noise2 = noise2_param*randn(N,N);
if noise2_param~=0
    SNR2 = 20*log10(norm(X(:),'fro')./norm(noise2(:),'fro'))
end
[u, s, v] = svd(X+noise2);
tsvd_rank = numel(find(abs(diag(s)) >= tsvd_thr));
compression_ratio = tsvd_rank*(1+2*N)./(N*N);
fprintf('Approximation has rank=%d \nCompression Ratio: %2.1f %%\n',...
         tsvd_rank,100*compression_ratio)
%T2 will be reactivated every 'T2_skip' iterations.
T2_skip = 5;
%Generate data:

p2 = tsvd(X+noise2,tsvd_thr,'soft');
p2_show = tsvd(X+noise2,tsvd_thr);

%T1 promotes sparsity in standard basis on first component,
%and sparsity in DCT basis in second component. 
dct_thr = 45;
img_thr = 10;
T1 = @(x) T1_subf(x,[img_thr, dct_thr]);
  
%Initialize variables for loop
xold = zeros(N,N,2);
xold2=xold;
xold3=xold;
err = Inf;
tol = 10^-6;
max_it = 1000;
E = zeros(max_it,1);
k=0;
T2 = @(x) p2 + x - F2(x);
T3 = @(x) sat(x,255,0);

fprintf('\n Starting main. alg...\n')  
k=0;
  
R1x = zeros(N,N,2,2);
R1x(:,:,:,1) = T1(xold);
R1x(:,:,:,2) = T2(xold);
nn = 1/2;
while k<max_it
    R1x(:,:,:,1) = T1(xold);
    if mod(k,T2_skip(1))==0 
    R1x(:,:,:,2) = T2(xold);
    end
    xnew = T3(nn*(sum(R1x,4)));

    xold=xnew;

    k=k+1;
    if mod(k,50)==0
      fprintf('%2.1f %%\n',100*k/max_it)
    end
end
  
figure(1)
subplot(2,2,1)
imshow(uint8(X))
title('Original')

subplot(2,2,2)
imshow(uint8(p2_show))

subplot(2,2,3)
imshow(uint8(xnew(:,:,1)))

subplot(2,2,4)
imshow(uint8(xnew(:,:,2)))

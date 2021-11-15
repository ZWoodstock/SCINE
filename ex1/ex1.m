close all; clear all;
pkg load image
page_output_immediately(1)

x_true = double(imread('cameraman.png'));

nxt = norm(x_true,'fro');

windowsize = 15; %Windowsize should be an ODD number
sigma = 3.5;% 
%This is the modelled saturation threshold.
sat_thr = 60;
ker = fspecial('gaussian',windowsize,sigma);
ker = ker./sum(sum(ker));
F1 = @(x) imfilter(isat(imfilter(x,ker,'conv','replicate'),...
          sat_thr),ker,'conv','replicate');
R1 = @(x) isat(imfilter(x,ker,'conv','replicate'), sat_thr);
%This is the real saturation threshold
sat_thr_r = sat_thr;
ker = fspecial('gaussian',windowsize,sigma);
F1_r = @(x) imfilter(isat(imfilter(x,ker,'conv','replicate'),...
	sat_thr_r),ker,'conv','replicate');
R1_r = @(x) isat(imfilter(x,ker,'conv','replicate'),sat_thr_r);

%T1 will be reactivated every 'T1_skip' iterations.
T1_skip = [1];

%uniform noise:
noise1_param = 20;
noise1 = unifrnd(0,noise1_param,size(x_true));
x1_corrupted = sat(x_true+noise1,255,0);
noise1 = x_true - x1_corrupted;
if noise1_param~=0
    SNR1 = 20*log10(norm(x_true,'fro')./norm(noise1,'fro'))
end

noise2_param = 1;
noise2 = unifrnd(0,noise2_param,size(x_true));
x2_corrupted = sat(x_true+noise2,255,0);
noise2 = x_true - x2_corrupted;
if noise2_param~=0
    SNR2 = 20*log10(norm(x_true,'fro')./norm(noise2,'fro'))
    temp = angle(fft2(x_true+noise2));
else
    temp = angle(fft2(x_true));
end

data = struct('phase',temp,'h',exp(i*temp));
T2 = @(x) proj_phase(x,data);
T2_skip = [1];

%Generate data:

p1_show = scal_img(R1_r(x_true+noise1));
p1 = imfilter(scal_img(p1_show,0,255,[0,sat_thr]),...
     ker,'conv','replicate');

%Initialize variables for loop

xold = zeros(size(x_true));
xold2=xold;
xold3=xold;
err = Inf;
tol = 10^-6;
max_it = 1000000;
E = zeros(max_it,1);
k=0;
T1 = @(x) p1 + x - F1(x);
noise4_param = 2;
if noise4_param~=0
noise4 = unifrnd(0,noise4_param,size(x_true));
else
noise4 = zeros(size(x_true));
end
p4 = mean(mean(x_true+noise4)); 
T4 = @(x) x - mean(mean(x)) + p4;

T3 = @(x) sat(x,255,0);

n_obs = 3; nn = 1./n_obs;
% 
max_it1 = 100000;
fprintf('Checking for inconsistency of the problem...\n')
xoldp = xold;
B = 20*nxt;
INC_flag = 0;
while err>tol && k<max_it1
    T1x = T1(xoldp); n1 = norm(xoldp-T1x);
    T2x = T2(xoldp); n2 = norm(xoldp-T2x);
    T3x = T3(xoldp); n3 = norm(xoldp-T3x);
    T4x = T4(xoldp); n4 = norm(xoldp-T4x);
    d = nn*(T1x+T2x+T3x+T4x) - xoldp;
    %Extrapolation Step
      if mod(k,3)==2
        lambda = 0.5*(nn*(n1.^2 + n2.^2+ n3.^2+n4.^2 ))./(norm(d).^2);
      else
        lambda = 1.99*(nn*(n1.^2 + n2.^2+ n3.^2+n4.^2 ))./(norm(d).^2);
      end     
    xnewp = xoldp + lambda*d;
    err = norm(xnewp-xoldp,'fro');
    nxn= norm(xnewp,'fro');
    if  nxn>B
      fprintf('|x| blew up; inconsistency likely.\n')
      %Set error to zero to break out of loop
      err=0;
      INC_flag = 1;
    end
    if mod(k,500)==0
      fprintf('%d iterations;|x|/|x_true| =%2.1f %%; Residual=%1.2e\n',k,100*nxn/nxt,err)
    end
    xoldp=xnewp;
    k=k+1;
end


fprintf('\n Starting main. alg...\n')  
k=0;

[N M]=size(xold);
R1x = zeros(N,M,3);
R1x(:,:,1) = T1(xold);
R1x(:,:,2) = T2(xold);
R1x(:,:,3) = T4(xold);

nn2 = 1/3;
while k<max_it
    if mod(k,T1_skip(1))==0
      R1x(:,:,1) = T1(xold);
    end
    if mod(k,T2_skip(1))==0
      R1x(:,:,2) = T2(xold);
    end
    R1x(:,:,3) = T4(xold);
    xnew = T3(nn2.*(sum(R1x,3)));

    E(k+1) = norm(xnew-x_true,'fro');

    xold=xnew;

    k+=1;
end
  
figure(1)
subplot(2,2,1)
imshow(uint8(x_true))
title('Original')

subplot(2,2,2)
imshow(uint8(p1_show))
title('Blurred, then saturated')

subplot(2,2,4)
imshow(uint8(xnew))
title('Recovery')

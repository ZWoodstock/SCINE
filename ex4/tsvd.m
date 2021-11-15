function out = tsvd(I,thr,method)
%Computes the truncated SVD of the input image I by
%hard-thresholding or soft-thresholding its singular values at the
%threshold 'thr'.
[u, s, v] = svd(I);
if nargin==3 && method=='soft'
    s = prox_abs(s,thr);
else
   %if the method is not specified to be soft-thresholding, we will
   %hard-threshold by default.
    s = (s.*(s>=thr));
end
out = u*s*(v');

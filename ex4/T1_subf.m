function out = T1_subf(x,gam)
%This is a subfunction for the operator T1 in the main script
%ex4.m
%Input: x (collated grayscale images);
%       gam (parameter for thresholding sparsity). Could be one
%       index or two, if you want to use separate thresholding
%       parameters for the DCT domain vs normal domain.
%Output: out = collated images, which have been thresholded in the
%standard basis (first component) and the DCT basis (second
%component).
  if numel(gam)==1
     gam(2) = gam;
  end
  out = zeros(size(x));
  out(:,:,1) = prox_abs(x(:,:,1),gam(1));
  out(:,:,2) = idct2(prox_abs(dct2(x(:,:,2)),gam(2)));
end


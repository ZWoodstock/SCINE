function out=scal_img(I,l,h,target)
%Rescale an image I to have its values between [0,255].
%If only one input is received, this will rescale an image 
%using its minimum and maximum values.
%If lower (l) and upper (h) bounds are received, we truncate all
%values outside of [l,h] and then rescale. 
%If the optional entry 'target' is used, then we rescale the image
%values between l and h to the target values specified in target.
%If target is scalar, we rescale to [0,target]; if target is a
%vector, we scale to [min(target), max(target)]. 
if nargin==1
    m = min(min(I));
    M = max(max(I));
    a = 255./(M-m);
    b = -m*a;
    out = a*I + b;
elseif nargin<=3
    lh = [min(l,h), max(l,h)];
    a = 255./(lh(2)-lh(1));
    b = -lh(1)*a;
    out = a*sat(I,lh(2),lh(1)) + b;
else
     lh = [min(l,h), max(l,h)];
     if numel(target)==1
     	target = [0, target];
     else
         target = [min(target), max(target)];
     end
     a = (target(2)-target(1))/(lh(2)-lh(1));
     b = target(1)- a*lh(1);
     out = a*sat(I,lh(2),lh(1)) + b;
  end
end

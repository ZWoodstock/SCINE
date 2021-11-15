function out = isat(I,thr)
%Inputs an image I (with entries between [0,255]) and a threshold
%'thr' in [0,255] and saturates the image. Saturation only occurs
%at the higher light levels (i.e. the image input will be truncated
%to be between [0, thr].
b = (I<=thr);
out = I.*b + thr.*(~b);
end

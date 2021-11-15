function out = proj_phase(x,data)
%This function inputs a 2-D (potentially color) image x, and
%outputs its projection onto the cone of functions with Fourier
%phase equal to the phase prescribed in data in the specified
%bandwidth. The input "data" must be a structure with at least two
%inputs:
% 'data.phase' must be the angle of the phase of each component of
%            a Fourier transform (same size as x)
%  'data.h' must be exp(i*data.phase), precomputed.
%If data.channel is specified (for use with color images) then this
%function computes the projection onto the specified phase only in
%the channel specified, leaving other channels of x untouched.
if isfield(data,'channel')==0
    data.channel = 1;
end
x_fft = fft2(x(:,:,data.channel));
out = x;
out(:,:,data.channel) = real(ifft2(abs(x_fft).*max(0,cos(data.phase - angle(x_fft))).*data.h));
end

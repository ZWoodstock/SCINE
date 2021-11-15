function out=sat(x,sens_up,sens_dn)
%sat(x,sens_up,sens_dn)
%Truncates a 1-D signal x between thresholds sens_up and sens_dn.
%If only two inputs are entered, we assume the interval to be
%symmetric
if nargin==2
   sens_up = abs(sens_up);
   sens_dn = -sens_up;
end
out = sens_up.*(x>sens_up) + x.*(x<=sens_up).*(x>=sens_dn) +...
	sens_dn.*(x<sens_dn);
end

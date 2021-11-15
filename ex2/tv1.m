function out = tv1(x,adj)
%Performs the total variation operation on x; optionally performs
%the adjoint as well.
%tv(x) = (x1-x2, x2-x3, ..., xn-1 - xn)
%note that \|L\|= 2 and \|L\|^2 = 4. so rescale by 0.25 to make
%tv1(x,1) firmly nonexpansive.
if nargin<2
%compute tv(x)
    out = diff(-x(:));
elseif (adj==1 || strncmp(adj,'adj',3))
%do not perform adjoint unless asked to do so.
%computes tv_adjoint(tv(x))
    out = [out(1); diff(out);-out(end)];
elseif strncmp(adj,'only',4)
%computes tv_adjoint(x)
    out = [x(1); diff(x); -x(end)];       
end
end

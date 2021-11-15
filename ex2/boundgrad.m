function out = boundgrad(x,gam)
%Computes the firmly nonexpansive operator  0.25*L^*(Id - proj_C) L,
%where L is the total variation operator in 1d, and C is the
%infinity-ball at 0 of radius 'gam'.
out = tv1(x);
%proj_C(x) = sat(x,gam,-gam);
out = (0.25).*tv1(out - sat(out,gam,-gam),'only');
end


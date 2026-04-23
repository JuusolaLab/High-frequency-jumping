function [BumpShape BumpDuration] = Bumpcalc(A,n,tao,tt)
%BUMPCALC Summary of this function goes here
%   Detailed explanation goes here
BumpShape = (( tt./tao ).^n) .* exp(-tt./tao)/tao /gamma(n+1)*tt(2);
I_in                = BumpShape;
[BumpMax, BumpMaxIn]= max(I_in);
II                  = I_in(BumpMaxIn+1:end);
II_ind              = find(II<0.02);
BumpDuration        = II_ind(1) + BumpMaxIn;
BumpShape = A/trapz(tt,BumpShape)*BumpShape;
end


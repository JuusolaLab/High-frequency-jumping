function [BumpShape BumpDuration] = Bumpcalc(A,n,tao,tt)
%BUMPCALC Calculates the shape of gamma bump. 
%A area of bump
%n gamma shape of bump
%tao gamma time constant of bump
%tt time vector
BumpShape = (( tt./tao ).^n) .* exp(-tt./tao)/tao /gamma(n+1)*tt(2);
I_in                = BumpShape;
[BumpMax, BumpMaxIn]= max(I_in);
II                  = I_in(BumpMaxIn+1:end);
II_ind              = find(II<0.02);
BumpDuration        = II_ind(1) + BumpMaxIn;
BumpShape = A/trapz(tt,BumpShape)*BumpShape;
end


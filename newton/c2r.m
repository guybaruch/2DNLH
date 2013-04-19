function r = c2r(c)
% C2R  reshapes an NMx1 complex vector as an 2NMx1 real vector
    global T;
	r(1:2:T.NM2-1,1) = real(c);
	r(2:2:T.NM2,1) = imag(c);
end % c2r()

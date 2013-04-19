function c=r2c(r)
% R2C  reshapes an 2NMx1 real vector as an NMx1 complex vector
    global T;
	c = r(1:2:T.NM2-1,1) + r(2:2:T.NM2,1)*j;
end % r2c()

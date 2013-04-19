function transverse_Laplacian_3D_2()
% TRANSVERSE_LAPLACIAN_3D_2  Creates the second-order transverse Laplacian on a
% Cylindrical grid, for both Dirichlet BCs and local Sommerfeld BCs.
	
	global T;
	
	if (T.half_int_xs==0) error('radial only uses half_int_xs'); end;
		
	e = ones(T.M,1);
	rhos = ( (1:T.M)-0.5 )' ;
	l = ones(T.M,1)./rhos;

	% spdiags( [ a b c], -1:1, ...)
	% b1 c2
	% a1 b2 c3
	% -  a2 b3 c4
	%
	% so easier to define L', then invert
	dgs_t = e*[1 -2 1]+l*[1 0 -1]/2;
	
	% remove c1, a_M (outside "matrix")
	dgs_t(1,3)=0;
	dgs_t(T.M,1)=0;
	
	
	%T.L_x = spdiags(dgs, -1:1, T.M, T.M)';
	%full(T.L_x)

	% symmetric BC at axis : no correction
	% B0 = [0 0];
	
	if ( T.trans_SBC==0)
		% Dirichlet BC
		B1 = - T.M/(T.M-1/2);
	else
		% Sommerfeld BC
		alpha = -T.k0 * besselh(1,1,T.k0*T.X_max)...
				/ besselh(0,1,T.k0*T.X_max); 

		alpha_h = alpha*T.h_x;

		B1 = T.M/(T.M-1/2) * (2+alpha_h)/(2-alpha_h);
	end;
	
	dgs_t(T.M,2) = -2 + B1;


	T.L_x_dgs_t = dgs_t;
	T.L_x = spdiags(T.L_x_dgs_t, -1:1, T.M, T.M).';
	%full(T.L_x)
	T.L_x = T.L_x/T.eta_x^2;


end % transverse_Laplacian_3D_2();

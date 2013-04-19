function init_long()
% INIT_LONG  creates the longitudinal Laplacian operator.
% Assumes that init_long_grid and init_trans were called.

	global T;

	T.k_pars = sqrt( ones(T.M,1) - T.lambdas )*T.k0;
	T.eta_pars = T.k_pars * T.h_z;

	if (2==T.order)	% no anodes
		% prepare the line [1 -2 1] and the matrices G,F 
		longitudinal_Laplacian_2();

	elseif (4==T.order) 
		
		if (T.compact==0)
			% prepare the line [-1 16 -30 16 -1] and the matrices G_1,2,3,4 ,F_0,1 
			longitudinal_Laplacian_4();
		elseif (T.compact==1)
			% prepare the line [1 -2 1]/(1+etas^2/12) and the matrices G,F
			longitudinal_Laplacian_4_compact();
		end;

	end;

end % init_long()

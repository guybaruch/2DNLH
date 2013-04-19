function longitudinal_Laplacian_2()
% LONGITUDINAL_LAPLACIAN_2  Creates the second-order 1D Laplacian L_z
% 	Using the standard centered-difference [1 -2 1] stencil
%	q^l and the matrices G, F
%	the matrices B_0,1  and Psi^-1 F Psi.

	global T;
	
	T.L_z = [1 -2 1] / T.eta_z^2;

	app_direct = 1;		% 0: -2+q , 1: direct centered-difference 

	e = ones(T.M,1);

	if (app_direct)
		C = (e*2 + T.eta_pars*j) ./ (e*2 - T.eta_pars*j);
		T.G = diag(C);
		f = (-T.eta_pars*4*j) ./ (e*2 - T.eta_pars*j);
		T.F = diag(f);
	else
		t2 = T.eta_pars.*T.eta_pars;

		q = e - t2/2 + sqrt(e-t2/4).*T.eta_pars*j ;
		
		T.G = diag(q);

		sq = sqrt(q);
		T.F = diag(sq.*(conj(q)-q));
	end;

	T.Bz = U2E(T.G);
	T.Bf = U2E(T.F);
end

function longitudinal_Laplacian_4_compact()
% LONGITUDINAL_LAPLACIAN_4_COMPACT  Creates the fourth-order 1D Laplacian L_z
% 	in the semi-compact case.
	global T;

	if (T.compact==0) error(); end;

	T.L_z = [1 -2 1] *(T.eta_z^-2 + 1/12);

	e = ones(T.M,1);

	L1 = e*(1+(T.eta_z^2)/12);
	L0 = L1 - T.eta_pars.^2/2;
	r = L0./L1;
	q = r +j*sqrt(-r.^2+1);
	%q = r +sqrt(r.^2-1);
	
	T.G = diag(q);
	T.F = diag(q.^-T.interface_distance.*(e./q-q));

	T.Bz = U2E(T.G)*( 1+(T.eta_z^2)/12 );
	T.Bf = U2E(T.F)*( 1+(T.eta_z^2)/12 );


	%T.L_z_interface_line = [-1 8 -14 8 -1] / (6*T.eta_z);
	T.L_z_interface_line = [4 -27 108 -170 108 -27 4] / (66*T.eta_z);
	
end

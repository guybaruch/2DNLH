function init_trans()
% INIT_TRANS  Creates the transverse Laplacians T.L_x
%   and their diagonal elements T.L_x_diags transposed, 
%   i.e. T.L_x = spdiags(T.L_x_dgs_t).' .
%
% Assumes a previous call to init_trans_grid

% spdiags( [ a b c], -1:1, ...)
% b1 c2
% a1 b2 c3
% -  a2 b3 c4

% spdiags( [ a b c d e] , -2:2,...)
% c1 d2 e3
% b1 c2 d3 e4
% a1 b2 c3 d4 e5
% -  a2 b3 c4 d5 e6
	global T;


	if (T.radial && T.x_asymmetric)
		error('no radial and x_asymmetric.');
	end;
	if (T.order==2 && T.compact)
		error('no SO and compact.');
	end;
	
	switch (T.radial + T.order)
		case 2 % O(h^2) cartesian
			transverse_Laplacian_2D_2();
		case 3 % O(h^2) radial
			transverse_Laplacian_3D_2();
		case 4 % O(h^4) cartesian
			transverse_Laplacian_2D_4();
		case 5 % O(h^4) radial
			transverse_Laplacian_3D_4();
		otherwise
			error('init_trans: radial+order mismatch')
	end;

	% guybar : later recreate the save/restore mechanisms ?
	%[T.Psi, T.Lambda] = eigs(T.L_x, T.M);
	[T.Psi, T.Lambda] = eig(full(T.L_x));

	T.lambdas = -diag(T.Lambda);

	%[ T.Psi_L, T.Psi_U, T.Psi_P ] = lu(T.Psi);
	T.Psi_inv = inv(T.Psi);

end % init_trans()

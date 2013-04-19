function transverse_Laplacian_2D_2()
% TRANSVERSE_LAPLACIAN_2D_2  Creates the second-order transverse Laplacian on a
% Cartesian grid, for both Dirichlet BCs and local Sommerfeld BCs.
	global T;

	% spdiags( [ a b c], -1:1, ...)
	% b1 c2
	% a1 b2 c3
	% -  a2 b3 c4
	%
	% so easier to define the transposed diags , then spdiags and transpose

	e = ones(T.M,1);

	dgs_t = e*[1 -2 1];
	
	% remove c1, a_M (outside "matrix")
	dgs_t(1,3)=0;
	dgs_t(T.M,1)=0;


	if ( T.trans_SBC==0)
		% Dirichlet BC
		if (T.half_int_xs)
			B1 = -1;
		else
			B1 = 0;
		end;
	else
		q_x = exp(j*T.eta_x);
		%q_x = (2+T.eta_x*j)/(2-T.eta_x*j);
		%q_x = ( 2 - T.eta_x^2 + j*T.eta_x*sqrt(4-T.eta_x^2) ) / 2 ; 
		% Sommerfeld BC
		B1 = q_x;
	end;
	% T.L_x(T.M,T.M) = -2 + T.L_x_B1;
	dgs_t(T.M,2) = -2 + B1;

	
	% BC at lower x (0 or -X_max)
	if (T.x_asymmetric)
		% same BC in -X_max as in X_max
		B0 = [0 B1];
	else
		% symmetric BC at axis 
		if (T.half_int_xs)
			B0 = [1 0];
		else
			B0 = [0 1];
		end;
	end;
	
	% [ b1 a1]
	dgs_t(1,1:2) = dgs_t(1,1:2) + B0 ;
	% T.L_x(1,1:2) = [-2 1] + T.L_x_B0;

	T.L_x_dgs_t = dgs_t;
	T.L_x = spdiags(T.L_x_dgs_t, -1:1, T.M, T.M).';
	%full(T.L_x);
	T.L_x = T.L_x/T.eta_x^2;
	
end % transverse_Laplacian_2D_2();

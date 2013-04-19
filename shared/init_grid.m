% function init_grid()
%
% creates the x-z grids

function init_grid(k0, Z_max, N, X_max, M, x_asymmetric, half_int_xs, SBC, order, compact)
    
	global T;
	
	T.k0 = k0;
	T.order = order;
	T.compact = compact;

	init_trans_grid(X_max, M, x_asymmetric, half_int_xs, SBC);
	init_long_grid(Z_max, N);

	T.h_r = T.h_z/T.h_x;    % h ratio : L_x*hr2 +Lz

end % function init_grid()

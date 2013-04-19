function init_trans_grid(X_max, M, x_asymmetric, half_int_xs, SBC)
% INIT_TRANS_GRID  Creates the grid in the transverse direction.
	global T;
	
	T.x_asymmetric = x_asymmetric;
	T.half_int_xs = half_int_xs; 
	T.M = M;	
	T.X_max = X_max; 	
	T.trans_SBC = SBC;
	
	T.ms = 0:M-1;
	
	if (half_int_xs)
		if (x_asymmetric)
			% x_1=-X_max+dx/2, .. , x_M = X_max-dx/2
			T.h_x = 2*X_max/M;
			T.xs = (T.ms+(1-M)/2)*T.h_x;
		else
			% x_1=dx/2, .. , x_M = X_max-dx/2
			T.h_x = X_max/M;
			T.xs = (T.ms+1/2)*T.h_x;
			
		end;
	else
		if (T.trans_SBC || T.x_asymmetric)
			error('transverese integer nodes with SBC or x_asymmetric unimplemented');
		end;

		% x_1=0, .. , x_M = X_max-dx
		T.h_x = X_max/M;
		T.xs = (T.ms)*T.h_x;
	end;

	T.eta_x = T.h_x * T.k0;

end % init_trans_grid()

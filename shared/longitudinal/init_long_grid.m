function init_long_grid(Z_max, N)
% INIT_LONG_GRID  Creates the uniform grid (integer or semi-integer) in the 
% longotudinal direction.

	global T;

	T.Z_max = Z_max;
	T.h_z = Z_max/N;
	
	T.vN = N;

	%: half_int, vd=2: N=vN+2*vd = N+4
	% -3/2 -1/2 | 1/2 , ... , vN-1/2 | vN+1/2, vN+3/2
	% 1     2   | 3,    ... , vN+2,  | vN+3,   vN+4
	% 0     1   | 2,    ... , N-3,   | N-2,    N-1

	%: int, vd=4: N = vN+2vd+1 = vN+9
	% -4, ... , -1 | 0 | 1 , ... , vN-1 | vN  | vN+1, ..., vN+4
	% 1 , ... ,  4 | 5 | 6 , ... , vN+4 | vN+5| vN+6, ..., vN+9
	% 0 , ... ,  3 | 4 | 5 , ... , N-6  | N-5 | N-4,  ..., N-1

	if (2==T.order)	% no anodes
		T.vd = 0;
		T.half_int_zs = 1;
	elseif (4==T.order) % +2 anodes from each side
		if (T.compact==0)
			T.vd = 2;
			T.half_int_zs = 1;
		else
			T.vd = 4;
			T.half_int_zs = 0;
			T.interface_nodes = [4 T.vN+4];
			T.interface_distance = 4;  % how far is the interface from n=0 ?
		end;
	end;

	T.int_zs = 1-T.half_int_zs;

	T.N = N+2*T.vd + T.int_zs;
	T.vns = (0:T.N-1) ... % base list
			- T.vd ... % additional nodes
			+ T.half_int_zs/2 ; % +1 node for integer nodes

	T.zs = T.vns*T.h_z;

	T.n_int_start = 0+T.vd+T.half_int_zs;
	T.n_int_end = T.vN+T.vd;

	if (2==T.order)	% no anodes
		T.p_ind_int = 1:T.N*T.M;
	elseif (4==T.order) % +2 anodes from each side
		if (T.compact==0)
			T.p_ind_int = 1+T.M*2:(T.N-2)*T.M;
		else
			T.p_ind_int = 1+T.M*5:(T.N-6)*T.M;
		end;
	end;

	T.eta_z = T.h_z * T.k0;

	T.NM = T.M*T.N;
	
end % init_long_grid()

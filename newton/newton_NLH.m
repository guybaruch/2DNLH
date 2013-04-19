% function newton_NLH
%
% solves tha NLH using the newton iterations.
%	(real) K0, X_rp, Z_rp, xmax, ymax, p=N/Nc,
%	(int) max_I, max_J, flags=FR4S
				

function newton_NLH(radial, sigma, x_asymmetric, compact, k0, Z_max, N, X_max, X_shift, M,...
					epsilon, max_I, IG, d_tau, IC)
% solves 2D the NLH using the newton iterations.

	disp(['newton_NLH(' ...
		  'radial=' int2str(radial) ...
		  ', sigma=' int2str(sigma) ...
		  ', x_asymmetric=' int2str(x_asymmetric) ...
		  ', compact=' int2str(compact) ...
		  ', k0=' num2str(k0) ...
		  ', Z_max=' num2str(Z_max) ...
		  ', N=' int2str(N) ...
		  ', X_max=' num2str(X_max) ...
		  ', X_shift=' num2str(X_shift) ...
		  ', M=' int2str(M) ...
		  ', epsilon=' num2str(epsilon) ...
		  ', max_I=' int2str(max_I) ...
		  ', IG=' int2str(IG) ...
		  ', d_tau=' num2str(d_tau) ...
		  ', IC="' IC   ...  
		  '")']);
	
	global T;


	tic
	addpath('../shared', '../shared/longitudinal', '../shared/transverse');
	
	[SBC order] = newton_NLH_paprams(radial, sigma, epsilon, X_shift, max_I, IG, d_tau,IC);
	
	format compact;

	format short e;

	init_grid(k0, Z_max, N, X_max, M, x_asymmetric, T.half_int_xs, SBC, order, compact);
	
	generate_E_inc(k0, IC);

	build_L();


	prep_test();

	get_initial_E();
	
	build_J0(); T.J=T.J0;
	T.dV2 = zeros(T.NM2,1);
	
	T.V = reshape(T.E, T.NM, 1);
	disp('Building time: '); 	
	toc
	tic
	for iter = 1:max_I
		
		newton_NLH_step(iter);
		
		if (T.max_dV > 1000)
			display('Divergence likely, max |dE|>1000')
			break
		end;

		dE_thresh = 1e-2;
		if (T.max_dV < dE_thresh && T.d_tau ~= 1)
			display(['convergence likely, max |dE|<' ...
					 num2str(dE_thresh) ', changing to d_tau=1']);
			T.d_tau = 1;
		end;
		
		if (T.max_dV < 1e-12 )
			display('convergence likely, max |dE|<1e-12, stopping');
			T.max_I = iter;
			break;
		end;


	end; % for iter
	
	toc
    
	if (max_I) 
		%		T.df(T.t_i,4:5) 
	end;

	if (T.trace_jumps==0)
		T.res = reshape(T.res_V, T.M, T.N);
		T.dE = reshape(T.dV, T.M, T.N);
	end;
		
	do_plots();
	do_saves();
	

end % newton_NLH

function [SBC order] = newton_NLH_paprams(radial, sigma, epsilon, X_shift, ...
										  max_I, IG, d_tau, IC)

	global T;
	
	T.radial = radial;
	T.epsilon = epsilon;
	T.sigma = sigma;
   	T.d_tau = d_tau;
	T.max_I = max_I;
	

	
	T.echo_args = 1;

	T.save = 1; % 0 no, 1 final, 2 every trace
	T.plot = 2; % 0 no, 1 final, 2 every trace
	T.trace_jumps = 3; % 0: no trace
	
	T.interface_E_inc_fix = 1;
	T.IC = IC;
	T.IC_shift = X_shift;
	
	switch (IG)
		case 0
			T.initial_guess = 'zero'; 
		case 1
			T.initial_guess = 'load';
		case 2
			T.initial_guess = 'keep';
		case 3
			T.initial_guess = 'NLS';
		otherwise
			error('IC!=[0123');
	end;

	
	T.half_int_xs = 1;
	SBC = 1;
	order = 4;

	if(1)
		T.LA_method = 'LU';
	else
		T.LA_method = 'GJ';
		T.LA_method_w = 1;
		%T.LA_method_w = -0.001;
		T.LA_method_eps = 1e-4;
	end;
end


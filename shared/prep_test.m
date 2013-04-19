%  function prep_test()
%
% --------------
% defines the sizes needed for a 1D test.
%
% parameters (in T):
%	method:
%		order=2/4,
%		compact= 
%			0: no compact
%			1: compact z
% 	Grid: 
%		z_max, Nz, int_z_nodes=0/1
% 	Phys: 
%		k0, mu, epsilon, u_inc
%	ploting:
%		plotit, N_plots, flags ?
%
% defines (in T):
% 	grid: 
%		Zs, dz, eta, eta2 
% 	phys sizes: 
%		lambda0, lambda1, epsilon1
% 	fields: 
%		u, du, u_prev, res, I_CM, u_CM
%	plot_arrs:
%		df
%
%-------------------------------------------------

function prep_test()
	
	global T;


	% fields
	%	T.u = zeros(T.Nz,1);    % solution
	%	T.u_prev = zeros(T.Nz,1);    % solution
	%	T.du = zeros(T.Nz,1);
	%	T.dau = zeros(T.Nz,1);
	%	T.I_CM=zeros(T.N_int,1);
	%	T.u_CM=zeros(T.N_int,1);
	%	T.res=zeros(T.Nz,1);
	%	T.f0 = zeros(T.Nz,1);
	T.save_double = 1;
	
	T.lambda0 = 2*pi/T.k0;
	T.lambda1 = T.lambda0 / sqrt(1+T.epsilon);
	
	if (T.echo_args)
		%k0 = T.k0

		date
		
		disp(sprintf('radial=%d    sigma=%d    SBC=%d   order=%d	compact=%d    half_int_xs=%d',  T.radial, T.sigma, T.trans_SBC, T.order, T.compact, T.half_int_xs));

		disp(sprintf('domain\t(Z,X)=(%g, %g)\ngrid\t(N,M)=(%g, %g)', ...
					 T.Z_max, T.X_max, T.N, T.M));
		disp( sprintf('resolution\t(eta_z, eta_x) = k0*(h_z, h_x) = \n\t%g*(%g, %g) = (%g, %g)',...
				 T.k0, T.h_z, T.h_x, T.eta_z, T.eta_x));
				
		disp(sprintf('resolving power\t(Z_rp,X_rp)=(%0.2g, %0.2g)', ...
					 T.lambda0 / T.h_z, T.lambda0 / T.h_x));

	
		%epsilon = T.epsilon

		%d_tau = T.d_tau 
		
		if (T.radial==0 && T.sigma==2)
			p = sqrt(2*T.epsilon/(3*pi))*T.k0
		end;
		if (T.radial==1 && T.sigma==1)
			p = (T.epsilon*T.k0^2)/7.4492
		end;
		
		disp(sprintf('lambda0 = %0.4g\tlambda1 = %0.4g', T.lambda0, T.lambda1));		
		%save_double = T.save_double
		
		T.normalize_grid = 2;  % 0 none 1 lambda_0 2 L_DF=2k_0r_0
		switch T.normalize_grid
			case 0
		    disp('grid not rescaled');
		    T.gs = 1;
			case 1
			disp('grid rescaled by Lambda_0');
			T.gs = T.lambda0;
			case 2
			disp('grid rescaled by L_DF = 2k_0');
            T.gs = 2*T.k0;
		end;
	
	% freezing 
	T.E0 =  zeros(T.M,T.N); % initial guess
	
    T.A = T.L;


	% Newton stuff

	T.NM2 = T.NM*2;

	T.I_2 = [1 0 ; 0 1 ];
	T.i_sigma2 = [ 0 1 ; -1 0 ];
	T.sigma3 = [1 0 ; 0 -1];

	T.J0 = sparse(T.NM2); T.J = sparse(T.NM2);
	
	T.f_2 = zeros(T.NM2,1);
	T.V2 = zeros(T.NM2,1);  % solution split to real and imaginary
	T.dV2 = zeros(T.NM2,1);
	T.V2_prev = zeros(T.NM2,1);
%	T.v = zeros(T.Nz2,1);





	% plots
	T.t_i=0;	

	if (T.trace_jumps) 
		T.N_plots=T.max_I/T.trace_jumps;
	else
		T.N_plots=1;
	end;
	if T.plot
		T.df=zeros(T.N_plots,4);
	end;
	
end

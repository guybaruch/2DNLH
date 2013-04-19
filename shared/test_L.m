% function test_L()
%
% tests the linear-helmholtz operator L created by build_L
%
% M=32; test_L(20,2,4*M,1,M,0,1,1,0,4,0)


function test_L(k0, Z_max, N, X_max, M, x_asymmetric, radial, half_int_xs, SBC, order, compact)

	global T;
	
	addpath('longitudinal', 'transverse');
	
	T.radial = radial;


	init_grid(k0, Z_max, N, X_max, M, x_asymmetric, half_int_xs, SBC, order, compact)
	
	loc_generate_E_inc();

	build_L();
	
	loc_generate_E_continuous();
    
	format compact
	x_asymmetric 
	radial
	compact
	half_int_xs
	SBC
	order
    
	T.res = (T.L * T.V_cont - T.f0)*T.h_z;
	T.res_E = reshape(T.res,T.M,T.N);
	res_max = max(abs(T.res))

	test_error = 1;
	plotit=1;

	solution_method = 2; % 0 - direct; 1: [lupq], 2: [lupqr], -1: none

	if (test_error) % test error, not just residual
		switch (solution_method)
			case 0
				T.V = T.L\T.f0;
			case 1
				[l,u,p,q] = lu(T.L);
				T.V = q*(u\(l\(p*T.f0)));
			case 2
				[l,u,p,q,r] = lu(T.L);
				T.V = q*( u\(l\( p*(r\T.f0) )) );
			case -1
				disp('no sol');
				T.V = -T.f0;
		end
        T.E = reshape(T.V, T.M, T.N);

		T.dV = T.V-T.V_cont;
        T.dE = reshape(T.dV, T.M, T.N);
		dV_max = max(abs(T.dV))
		% also: plot
		
		if (plotit)
			if(1)
				figure(1) ; surf(T.zs, T.xs, abs(T.E )); title('|E|'); shading flat;
				figure(2) ; surf(T.zs, T.xs, abs(T.E_cont )); title('|E_{cont}|');shading flat;
				figure(3); surf(T.zs, T.xs, abs(T.dE )); title('|dE|'); shading flat;
				figure(4); surf(T.zs, T.xs, abs(T.E)-abs(T.E_cont)); title('d|E|'); shading flat;
				%figure(5); surf(T.zs, T.xs, abs(T.res_E )); title('|res_E|'); shading flat;
			end;
			if (0) % Einc0
				   %				figure(10); plot(T.xs, abs(T.Einc0)); title('Einc0');
				   le = T.L_x * T.Einc0 ;
				   dle = le - T.Einc0*(le(4))/T.Einc0(4);
				   figure(11); semilogy(T.xs, abs( dle)); title('(L_x-\lambda)E_{inc}^0');
				   abs_dle = max(abs(dle()))
				   figure(12); plot(T.xs, abs( le )); title('L_x*E_{inc}^0');
				   de = diff(T.Einc0)/T.h_x;
				   me = (T.Einc0(1:T.M-1)+T.Einc0(2:T.M))/2;
				   %				figure(13); plot(T.xs(1:T.M-1)+T.h_x/2, imag( de .*conj(me) )  ); title('imag(\psi_x \psi^*)');
				   ae = abs( de - j*T.k0*me )';
				   %ggg = [ T.xs(1:T.M-1)+T.h_x/2 ; ae ];
				   %format long; double(ggg(:, (-5:-1)+T.M))
				   %format 
				   P = polyfit( T.xs(T.M-5:T.M-1)+T.h_x/2 -1, ae(T.M-5:T.M-1),3);
				   P(4)
				   
				   figure(14); plot(T.xs(1:T.M-1)+T.h_x/2, abs( de - j*T.k0*me )  ); title('\psi_x -ik_0 \psi');
			end; % if plot E_inc0
		end;
	end;
    
    
function loc_generate_E_inc()
	T.Einc0 = zeros(M,1);
	T.Einc1 = zeros(M,1);
	
	if (radial)
		if (1~=X_max) error('unimplemented'); end;
		if (20~=k0) error('unimplemented'); end;
		if (SBC)
			op = optimset('TolFun',1e-12);
			alpha = -k0 * besselh(1,1,k0*X_max) / ...
					besselh(0,1, k0*X_max);
			x  = fsolve(@(x) alpha*besselj(0,x)+x*besselj(1,x) , 7.85, op);

			%x  = fsolve('(x-20)/(x+20)-exp(-2*j*x)', 7.8/X_max, op);
			T.k_perp = x/X_max;
		else
			op = optimset('TolFun',1e-14);
			x  = fsolve('besselj(0,x)', 7.85, op);
			T.k_perp = x/X_max;
		end;
		T.Einc0 = besselj(0, T.xs * T.k_perp).';
	else
		if (SBC)
			if (1~=X_max) error('unimplemented'); end;
			if (20~=k0) error('unimplemented'); end;

			op = optimset('TolFun',1e-14);
			x  = fsolve('tan(x)+20*j/x', 7.85, op);
			%x  = fsolve('(x-20)/(x+20)-exp(-2*j*x)', 7.8/X_max, op);
			T.k_perp = x/X_max;
		else
			%T.k_perp = 1.5*pi/X_max;
			T.k_perp = 2.5*pi/X_max;
		end;
		
		T.Einc0 = cos(T.xs * T.k_perp).';
	end;
	k_perp = T.k_perp
end;

function loc_generate_E_continuous()
	k_par = sqrt(k0^2 -T.k_perp^2)
		%if (1) 	k_par = conj( k_par ); end; % ggg try
	
	psi = exp ( T.zs *j*  k_par );

        T.E_cont = T.Einc0 * psi;
        T.V_cont = reshape(T.E_cont, T.NM,1);
end;

end % test_L

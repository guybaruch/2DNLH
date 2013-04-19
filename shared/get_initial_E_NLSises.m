function get_initial_E_NLSises()
% constructs the initial guess from NLS or coupled NLS.

	global T;

	N = T.N;
	M = T.M;

	eta_z = T.eta_z;

	A = zeros(M,N);
	B = zeros(M,N);
	E_NLS = zeros(M,N);

	if (T.compact==0) error('non compact with NLS IG is not implemnted'); end;

	%int_nodes = (0:T.vN)+T.vd+1;

	zs = T.zs;
	qs = exp(zs*j*T.k0);
		
	I_plus_L_z = speye(M) + T.L_x*i*T.eta_z/4;
	I_minus_L_z = speye(M) - T.L_x*i*T.eta_z/4;

	[L,U] = lu(I_minus_L_z);

	for ix = 1:3
		coupled_NLS_pass()
	end;


	T.A = A;
	T.B = B;
	T.E_NLS_IG = E;
	T.E = E;

	function coupled_NLS_pass()
			
		% A
		mx = 0;
		Psi_prev = T.Einc0;
		Psi = T.Einc0;
		for n=1:N
			nB = N-n;
			if (nB)
				anti_Psi = abs2( (B(:,nB) + B(:,nB+1))/2 );
			else
				anti_Psi = abs2( B(:,nB+1) );
			end;	
			step(n);
			mx = max(mx,max(abs2(A(:,n)-Psi)));
			A(:,n) = Psi;
		end;
		
		%disp(['max A=' num2str(mx)])

		% B
		Psi_prev = T.Einc1;
		Psi = T.Einc1;
		for n=1:N
			nA = N-n;
			if (nA)
				anti_Psi = abs2( (A(:,nA) + A(:,nA+1))/2 );
			else
				anti_Psi = abs2( A(:,nA+1) );
			end;	
			step(n);
			mx = max(mx,max(abs2(B(:,n)-Psi)));
			B(:,n) = Psi;
		end;
		disp(['max A,B=' num2str(mx)])

		for m=1:M
			E(m,:) = A(m,:).*qs + B(m,N:-1:1).*qs(N:-1:1);
		end
		
		
		function step(n)
			Psi_nph_1 = (Psi*3-Psi_prev)/2;
	
			Psi_prev = Psi;
	
			RHS1 = I_plus_L_z * Psi;
			
			RHS = RHS1+(i*T.epsilon*eta_z/2)*Kerr(Psi_nph_1,anti_Psi,n);
			Psi_pred = U\(L\RHS);
	
			Psi_nph_2 = (Psi+Psi_pred)/2;
			RHS = RHS1+(i*T.epsilon*eta_z/2)*Kerr(Psi_nph_2,anti_Psi,n);
			Psi  = U\(L\RHS);
	
			function kerr = Kerr(Psi_nph,anti_Psi,n)
				I = Psi_nph.*conj(Psi_nph);
				kerr = (I+2*anti_Psi).*Psi_nph;
				% saturate
				kerr = kerr ./(1+0.3*T.epsilon*abs2(kerr));
				if (n<T.vd || n>T.N-T.vd) 
					kerr = kerr*0;
				end
			end
	
			T.Psi = Psi;
		end; % step
	

	end; % function coupled_NLS_pass()


	

end %get_initial_E;
	

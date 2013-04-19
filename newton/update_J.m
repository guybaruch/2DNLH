function update_J()
% UPDATE_J  Updates the Jacobian at each Newton iteration.
    global T;

	if T.compact 
		V = T.V;
	
		I = V.*conj(V);
		VV = V.*V;
		
		switch (T.sigma)
			case 1
				K_E = 2*I;   % \partial (Kerr) / \partial E 
				K_Ec = VV;  % \partial (Kerr) / \partial E*
			case 2
				K_E = 3*I.^2;
				K_Ec = 2*I.*VV;
			case 3
				K_E = 4*I.^3;
				K_Ec = 3*I.*I.*VV;
			otherwise
				error('T.sigma != 1,2,3');
		end;
		
		K_E = K_E*T.epsilon;
		dgs = zeros(T.NM2,1);
		dgs(1:2:T.NM2-1) = real(K_E);
		dgs(2:2:T.NM2) = real(K_E);
		diag_K_E = spdiags(dgs, 0 ,T.NM2,T.NM2);

		K_Ec = K_Ec*T.epsilon;
		dgs = zeros(T.NM2,3);
		dgs(1:2:T.NM2-1,2) = real(K_Ec);
		dgs(2:2:T.NM2,2) = -real(K_Ec);
		dgs(1:2:T.NM2-1,1) = imag(K_Ec);
		dgs(2:2:T.NM2,3) = imag(K_Ec);
		diag_K_Ec = spdiags(dgs, -1:1 ,T.NM2,T.NM2);

		%s3 = ones(T.NM2,1); s3(2:2:T.NM2)=-1;
		
		T.J = T.J0 + T.dJdI * ( diag_K_E + diag_K_Ec);
			
	else
	% assume V==r2c(V2)
	m = [0 0 ; 0 0];
	
	V = T.V;
	if (2==T.n_int_start )
		V(1:2*T.M) = 0;
		V( (1-2*T.M:0)+T.NM  ) = 0;
	end;
	
	I = V.*conj(V);
	VV = V.*V;
	
	switch (T.sigma)
		case 1
			K_E = 2*I;   % \partial (Kerr) / \partial E 
			K_Ec = VV;  % \partial (Kerr) / \partial E*
		case 2
			K_E = 3*I.^2;
			K_Ec = 2*I.*VV;
		case 3
			K_E = 4*I.^3;
			K_Ec = 3*I.*I.*VV;
		otherwise
			error('T.sigma != 1,2,3');
	end;
	
	sigma3 = [1 0 ; 0 -1];
	epsilon = T.epsilon;
	
	if (0) % this is BAD performance wize
		for p=T.p_ind_int % go over all diagonal 2x2 blocks with Kerr
			ri = p*2-1; ii=p*2;
			K_Er = K_E(p); K_Ecr = K_Ec(p);
			m = [K_Er 0; 0 K_Er] + c2_2x2( K_Ecr )*sigma3;
			
			%T.J(ri:ii,ri:ii) = T.J0(ri:ii,ri:ii) + m*epsilon;
			m1 = T.J0_blkdiag(:,:,p) + m*epsilon;
			
			T.J(ri:ii,ri:ii) = m1;
		end; % for i
	else
			
		J_blkdiag = zeros(T.NM2+1,3);
		
		for p=T.p_ind_int 
			K_Er = K_E(p); K_Ecr = K_Ec(p);
			%m = [K_Er 0; 0 K_Er] + c2_2x2( K_Ecr )*sigma3;
			%m = [K_Er 0; 0 K_Er] + ...
			%	[real(K_Ecr) -imag(K_Ecr); imag(K_Ecr) real(K_Ecr) ]*sigma3;
			% [a -b ; b a]*sigma3 = [a b ; b -a]
			m = [K_Er+real(K_Ecr) imag(K_Ecr); imag(K_Ecr) K_Er-real(K_Ecr)];
			m=m*epsilon;
			J_blkdiag(p*2-1,2) = m(1,1);
			J_blkdiag(p*2,2) = m(2,2);
			J_blkdiag(p*2-1,1) = m(2,1);
			J_blkdiag(p*2,3) = m(1,2);
		end;
		T.J = T.J0 + spdiags( J_blkdiag, -1:1 ,T.NM2,T.NM2);
		
	end;
		
	end;

end % update_J


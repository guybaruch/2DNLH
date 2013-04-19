% function update_A
% given L and V, create A, which is the copy of L + the nonlinear terms

function update_A()
	global T;
	
	if ( (T.sigma~=1) && (T.sigma~=2) && (T.sigma~=3) ) error('T.sigma != 1,2,3'); end;

	I = T.V.*conj(T.V);
		
	if (2==T.sigma)
		I = I.*I;
	elseif (3==T.sigma)
		I = I.^3;
	end;
	
	I = I * T.epsilon;

	if (1)
		T.A = T.L + T.dAdI*spdiags(I, 0, T.NM, T.NM);
	else
		if (T.compact == 0)
			if (2==T.n_int_start )
				I(1:2*T.M) = 0;
				I( (1-2*T.M:0)+T.NM  ) = 0;
			end;
			
			T.A = T.L + spdiags(I, 0, T.NM, T.NM);
			
		else
			%I=I/12;
			V1 = I;
			ns = T.n_int_start;
			ne = T.n_int_end;
			V1(1:ns*T.M) = 0;
			V1(ne*T.M+1:T.NM) = 0;
			% [ n-1,m  n,m-1  n,m  n,m+1  n+1,m
			inds = [-T.M -1:1 T.M];
			dgs_t = ones(T.NM,1)*[1 0 -2 0 1]/12;
			% Kerr
			dgs_t(ns*T.M+1:ne*T.M,3)=dgs_t(ns*T.M+1:ne*T.M,3)+1;
			
			for n=ns:ne-1
				dgs_t(n*T.M+1:(n+1)*T.M,:) = ...
					dgs_t(n*T.M+1:(n+1)*T.M,:) + T.L_x_Laplacian_dgs_t*T.h_r^2;
	        end
			
			V5 = [ [zeros(T.M,1) ; V1(1+T.M:T.NM) ] ... % n-1,m
				   [ V1(2:T.NM) ; 0 ] ... % n,m-1
				   V1 ... % n,m
				   [0 ; V1(1:T.NM-1)] ... % n,m+1
				   [ V1(1:T.NM-T.M); zeros(T.M,1)] ];
			
			dgs_t_V = dgs_t.*V5;
			%figure(3); spy(dgs_t_V,'x')
			
			
			% add the interface nodes
			c1 = 3*T.eta_z/11;
	        dgs_t_V(1+(ns-1)*T.M:ns*T.M,3) = ...
	                dgs_t_V(1+(ns-1)*T.M:ns*T.M,3) + ...
	                I(1+(ns-1)*T.M:ns*T.M)*c1;
	        dgs_t_V(1+ne*T.M:(ne+1)*T.M,3) = ...
	                dgs_t_V(1+ne*T.M:(ne+1)*T.M,3) + ...
	                I(1+ne*T.M:(ne+1)*T.M)*c1;
	        
	
	        %T.s = spdiags(dgs_t_V, inds, T.NM, T.NM).';
			%figure(4); spy(s,'o')
			T.A = T.L + spdiags(dgs_t_V, inds, T.NM, T.NM).';
		end;
	end;
end

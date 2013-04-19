function build_J0()
% BUILD_J0  Builds the linear part of the NLH Jacobian
	% note: to see spdiags effect, do:
	%	A = spdiags( [1:10]'*[1:5], [-3 -1:1 3], 10,10);
	%	full(A)
	%	spdiags(A)

	global T;

	if (T.order~=4) error(''); end;

	if (T.compact)
		build_J0_compact();
	else
		build_J0_noncompact();
	end;

	function build_J0_noncompact()
		
		% spdiags( [ a b c d e f g h i j k ] , -5:5,...)
		% b <-> -4 ,d <-> -2, f<->0,  h<->2, j<->4
		%
		% f1 g2 | h3 i4 | j5 k6 |
		% e1 f2 | g3 h4 | i5 j6 |(k7) 
		% ------|-------|-------|-------|
		% d1 e2 | f3 g4 | h5 i6 | j7 k8 |
		% c1 d2 | e3 f4 | g5 h6 | i7 j8 |(k9)
		% ------|-------|-------|-------|-------|
		% b1 c2 | d3 e4 | f5 g6 | h7 i8 | j9 kA |
		% a1 b2 | c3 d4 | e5 f6 | g7 h8 | i9 jA |(kB)
		% ------|-------|-------|-------|-------|-------|
		% - (a2)| b3 c4 | d5 e6 | f7 g8 | h9 iA | jB kC |
		% -  -  | a3 b4 | c5 d6 | e7 f8 | g9 hA | iB jC | 
		% 
		% f goes to dgs(:,8)
		
		e = ones(T.NM2,1);
		%            -7->dgs(1)    0->dgs(8)    
		l_indices = [-4*T.M -2*T.M -5:5 2*T.M 4*T.M ];
		dgs = e * ( [-1 16 zeros(1,11) 16 -1]/(12*T.eta_z^2) );
		
		% these will be in the BC matrices
		dgs(1:2*T.M,2)=0;
		dgs(1:4*T.M,14)=0;
		dgs(T.NM2-T.M*4+1:T.NM2,2)=0;
		dgs(T.NM2-T.M*2+1:T.NM2,14)=0;
		
	
		% needed for update_J0
	    for p = 1:T.NM
			l = T.L(p,p);
			T.J0_blkdiag(:,:,p) = [real(l) -imag(l) ; imag(l) real(l)];
		end;
	
		
	    bwL = 2; % L_x bandwidth -bwL:bwL
		offs_L=1+bwL;
	    bwJ = 2*bwL+1; % J L_x bandwidth -bwJ:bwJ
		offs_J=1+bwJ;
		L_dgs = zeros(T.NM,offs_L+bwL);
		J_dgs = zeros(T.NM2,offs_J+bwJ);
		
		for k=-bwL:bwL
			L_dgs(1:T.NM-abs(k),offs_L+k)=diag(T.L,k);
			% remove BC matrix
			L_dgs(1:2*T.M-abs(k),offs_L+k)=0;
			L_dgs(T.NM-2*T.M+1-abs(k):T.NM,offs_L+k)=0;
		end;
		
		for k=-bwL:bwL
	
			J_dgs(1:2:T.NM2,offs_J+2*k) = real( L_dgs(:,offs_L+k));
			J_dgs(2:2:T.NM2,offs_J+2*k) = real( L_dgs(:,offs_L+k));
			
			if (k<0)
				J_dgs(2:2:T.NM2,offs_J+2*k+1) = -imag( L_dgs(:,offs_L+k));
				J_dgs(1:2:T.NM2,offs_J+2*k-1) = imag( L_dgs(:,offs_L+k));
			elseif (k==0)
				J_dgs(1:2:T.NM2,offs_J+2*k+1) = -imag( L_dgs(:,offs_L+k));
				J_dgs(1:2:T.NM2,offs_J+2*k-1) = imag( L_dgs(:,offs_L+k));
			else
				J_dgs(1:2:T.NM2,offs_J+2*k+1) = -imag( L_dgs(:,offs_L+k));
				J_dgs(2:2:T.NM2,offs_J+2*k-1) = imag( L_dgs(:,offs_L+k));
			end;
		end;
		
		for k=0:bwJ
			dgs(k+1:T.NM2,8+k) = J_dgs(1:T.NM2-k, offs_J+k);
			dgs(1:T.NM2-k,8-k) = J_dgs(1:T.NM2-k, offs_J-k);
		end;
		
		J = spdiags( dgs, l_indices, T.NM2, T.NM2 );
		
		B = build_J0_B(2);
	
		T.J0 = B+J;
	end; % function build_J0_noncompact
	

	function build_J0_compact()
	
		% first, the inner nodes
		l_indices = [-2*T.M -5:5 2*T.M];
	
		%
		% L: 1     2     3      4      5       6       7
		% J: (1)r (3)ir (5)ir (7)iri (9)ri (11)ri (13)r
		J0_dgs = zeros(T.NM2,13);
		
		Ldgs2Jdgs_real(1,1);
		Ldgs2Jdgs_real(7,13);
		for n=2:6
			Ldgs2Jdgs_complex(n,n*2-1);
		end;

        function Ldgs2Jdgs_real(L_diag,J_diag)
			J0_dgs(1:2:T.NM2-1,J_diag) = real(T.L_dgs(:,L_diag));
			J0_dgs(2:2:T.NM2,J_diag) = real(T.L_dgs(:,L_diag));
        end; % translate_real
		
        function Ldgs2Jdgs_complex(L_diag,J_diag)
        % for explanation, do spdiags( [1:4]'*[5:8] )
			Ldgs2Jdgs_real(L_diag,J_diag);
			J0_dgs(1:2:T.NM2-1,J_diag-1) = imag(T.L_dgs(:,L_diag));
			J0_dgs(2:2:T.NM2,J_diag+1) = -imag(T.L_dgs(:,L_diag));
        end; % translate_complex
		
		
		J = spdiags( J0_dgs, l_indices, T.NM2, T.NM2 );

		%figure(11); spy(J,'or'); title('J internal')
		%next, the interface nodes

		l_iface_inds = (-3:3)*2*T.M;
		dgs_iface_t = ones(8*T.M,1) * T.L_z_interface_line;
		L_interface_loc = spdiags(dgs_iface_t, l_iface_inds, 14*T.M, 8*T.M).';
		L_interface_loc = L_interface_loc(1+6*T.M:8*T.M,:);
					
		s1 = (T.interface_distance)*2*T.M;
		s2 = T.NM - 2*(s1+8*T.M);
		J_interface = blkdiag(...
			sparse(8*T.M,2*T.M), ...
			L_interface_loc, ...
			sparse((T.N-10)*2*T.M,(T.N-16)*2*T.M), ...
			L_interface_loc, ...
			sparse(8*T.M,2*T.M) );
		%figure(12); spy(J_interface,'ok'); title('J_interface')
		B = J_interface+build_J0_B(1);

		% this is so as not to add the J terms twice
		B(1:2*T.M,1:2*T.M) = B(1:2*T.M,1:2*T.M)-J(1:2*T.M,1:2*T.M);
		t1 = (T.N-1)*2*T.M;
		B(1+t1:T.NM2,1+t1:T.NM2) = B(1+t1:T.NM2,1+t1:T.NM2)...
			-J(1+t1:T.NM2,1+t1:T.NM2);

		
		T.J0 = B+J;
		%figure(13); spy(T.J0,'o'); title('J0')


		%now, build the dJdI
		dJdI_dgs = zeros(T.NM2,13);
		indices = [ -2*T.M -5:5 2*T.M  ];
		%dA2dJ(1,2);
			dJdI_dgs(1:2:T.NM2-1,1) = real(T.dAdI_dgs(:,1));
			dJdI_dgs(2:2:T.NM2,1) = real(T.dAdI_dgs(:,1));
		dA2dJ(2,3);
		dA2dJ(3,5);
		dA2dJ(4,7);
		dA2dJ(5,9);
		dA2dJ(6,11);
		%dA2dJ(5,9);
			dJdI_dgs(1:2:T.NM2-1,13) = real(T.dAdI_dgs(:,7));
			dJdI_dgs(2:2:T.NM2,13) = real(T.dAdI_dgs(:,7));
		
        function dA2dJ(ai,ji)
			dJdI_dgs(1:2:T.NM2-1,ji) = real(T.dAdI_dgs(:,ai));
			dJdI_dgs(2:2:T.NM2,ji) = real(T.dAdI_dgs(:,ai));
			dJdI_dgs(1:2:T.NM2-1,ji-1) = imag(T.dAdI_dgs(:,ai));
			dJdI_dgs(2:2:T.NM2,ji+1) = -imag(T.dAdI_dgs(:,ai));
		end;

		T.dJdI = spdiags(dJdI_dgs, indices, T.NM2, T.NM2);
		
		%figure(15); spy(dJdI_dgs,'ob'); title('dJdI\_dgs')
		
		%figure(16); spy(T.dJdI,'om'); title('dJdI')
		
	end; % function build_J0_compact

	function B = build_J0_B(w)
		% allocate the BC matrix
		
		rp = (T.N-w)*T.M;
		
		B0 = ones(T.M*2*w);
		B1 = ones(T.M*2*w);
		% fill BC matrix
		for i=1:T.M*w
			for j=1:T.M*w
				
				% lower zs
				l = T.L(i,j);
				C = [real(l) -imag(l) ; imag(l) real(l)];
				
				B0(i*2-1:i*2,j*2-1:j*2) = C;
	
				i1 = i+rp;
				j1 = j+rp;
				% upper  zs
				l = T.L(i1,j1);
				C = [real(l) -imag(l) ; imag(l) real(l)];
				
				B1(i*2-1:i*2,j*2-1:j*2) = C;
			end; % for j
		end; % for i

		
		t = (T.N-2*w)*T.M; 
		B = blkdiag(B0, sparse(t*2,t*2), B1);
	end;

end % build_J0 


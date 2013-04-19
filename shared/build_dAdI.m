function build_dAdI()
% BUILD_DADI  Creates the operator dAdI which is the nonlinear part of the
% implicit functional representation of the NLH.
% Specifically, if $F(E) = A(E)E - b$ and $L$ is the Laplacian then $L+dAdI$.
% See also: update_A
    global T;

	if (T.compact)
		if (T.order~=4) error(); end;
		%I=I/12;
		ns = T.n_int_start;
		ne = T.n_int_end;
		
		% [ n-1,m  n,m-1  n,m  n,m+1  n+1,m
		inds = [-T.M -2:2 T.M];
		dgs_t = ones(T.NM,1)*[1 0 0 -2 0 0 1]/12;

		% Kerr
		dgs_t(ns*T.M+1:ne*T.M,4)=dgs_t(ns*T.M+1:ne*T.M,4)+1;
		
		% remove outside and interface nodes
		dgs_t(1:ns*T.M,:) = 0;
		dgs_t(1+ne*T.M:T.NM,:)=0;

		%% L_x of Laplacian
		%for n=ns:ne-1
		%	dgs_t(n*T.M+1:(n+1)*T.M,:) = ...
		%		dgs_t(n*T.M+1:(n+1)*T.M,:) + ...
		%		[ zeros(T.M,1) T.L_x_Laplacian_dgs_t*T.h_r^2 zeros(T.M,1)];
        %end
		%%%%%% this .' save 98.5% of the cost ...
		% L_x of Laplacian
        dgs_t_L_x = ...
            [ zeros(T.M,1) T.L_x_Laplacian_dgs_t*T.h_r^2 zeros(T.M,1)].';
        tmp = dgs_t.';
        for n=ns:ne-1
            tmp(:,n*T.M+1:(n+1)*T.M) = tmp(:,n*T.M+1:(n+1)*T.M) + dgs_t_L_x;
        end;
        dgs_t = tmp.';

		% add the interface nodes
		c1 = 3*T.eta_z/11;
        dgs_t(1+(ns-1)*T.M:ns*T.M,4) = ...
        	dgs_t(1+(ns-1)*T.M:ns*T.M,4) + c1;
        dgs_t(1+ne*T.M:(ne+1)*T.M,4) = ...
        	dgs_t(1+ne*T.M:(ne+1)*T.M,4) + c1;
        
		T.dAdI_dgs_t = dgs_t;

        T.dAdI = spdiags(dgs_t, inds, T.NM, T.NM).';
		T.dAdI_dgs = spdiags(T.dAdI, inds);
		
%		figure(4); spy(T.NL,'o')
		
	else
		e = ones(T.NM,1);
		if (T.order==4)
			e(1:2*T.M)=0;
			e((T.N-2)*T.M+1:T.NM)=0;
		end;
		
		T.dAdI_dgs_t = e;

		T.dAdI = spdiags(e, 0, T.NM, T.NM);
	end;
end % build_dAdI

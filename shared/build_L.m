function build_L()
% BUILD_L  Creates the linear-Helmholtz part L of F_nm.
% See the definitions of F in [jcp 2009]
    global T;

    init_trans();
    init_long();

    hr2 = T.h_r^2;

    e = ones(T.NM,1);
    z = zeros(T.NM,1);

    T.L = sparse(T.NM);
    T.f0 = zeros(T.NM,1);
    
    if (2==T.order)
        dgs_t =  e * [ 1 0 -2 0 1];
        for n=0:T.N-1
            dgs_t(n*T.M+1:(n+1)*T.M,2:4) = ...
                dgs_t(n*T.M+1:(n+1)*T.M,2:4) + T.L_x_dgs_t*hr2;
        end;
        
        % add the linear-helmholtz
        dgs_t(:,3) = dgs_t(:,3) + T.eta_z^2;

        L = spdiags( dgs_t, [-T.M -1 0 1 T.M], T.NM, T.NM ).';
        
        % create the BC matrix
        z1 = zeros(T.NM,1); z1(1:T.M)=1; z1(T.NM-T.M+1:T.NM)=1;
        B = spdiags([z1], [0], T.NM, T.NM );
        
        % add the TWBC MxM blocks
        B(1:T.M,1:T.M) = T.Bz;
        B(T.NM-T.M+1:T.NM, T.NM-T.M+1:T.NM) = T.Bz;

        T.L = L+B;

        T.L = T.L/T.eta_z^2;

    else
        if (T.order~=4) error(); end;

        coef = 12*T.eta_z^2;
        
        if (T.compact==0)
    
            dgs_t = e * [-1 16 0 0 -30 0 0 16 -1];
            l_indices = [-2*T.M -T.M -2:2 T.M 2*T.M ];
    
            for n=0:T.N-1
                dgs_t(n*T.M+1:(n+1)*T.M,3:7) = ...
                    dgs_t(n*T.M+1:(n+1)*T.M,3:7) + T.L_x_dgs_t*hr2;
            end;
    
            % add the linear-helmholtz
            dgs_t(:,5) = dgs_t(:,5) + coef;
    
            L = spdiags( dgs_t, l_indices, T.NM, T.NM ).';
    
    
            % create the BC matrix
            t = (T.N-4)*T.M; 
            B = blkdiag(T.Bz0, sparse(t,t),T.Bz1);
        else
            % multiply by 12 b/c z is 1,-2,1 unlike x
            dgs_t = e * ( [1 0 0 -2 0 0 1] * (12+T.eta_z^2));
            l_indices = [-T.M -2:2 T.M];

            % add the linear-helmholtz
            dgs_t(:,4) = dgs_t(:,4) + coef;
            
            for n=0:T.N-1
                if (n==T.interface_nodes(1) || ...
                    n==T.interface_nodes(2))
                    % for the interface nodes, the previous are invalid
                    dgs_t(n*T.M+1:(n+1)*T.M,:) = ...
                        [ zeros(T.M,1) T.L_x_dgs_interface_t*coef zeros(T.M,1) ];
                else
                    dgs_t(n*T.M+1:(n+1)*T.M,2:6) = ...
                        dgs_t(n*T.M+1:(n+1)*T.M,2:6) + T.L_x_dgs_t*hr2;
                end;
            end;
    
            T.L_dgs_t = dgs_t;
            T.L_l_indices = l_indices;
            L = spdiags( dgs_t, l_indices, T.NM, T.NM ).';
            [tdgs, d] = spdiags(L);
            T.L_dgs = tdgs/coef;
            if (d'~=l_indices) error(''); end;
            
            % create the BC matrix
            t = (T.N-2)*T.M;
            B = blkdiag(T.Bz*12, sparse(t,t),T.Bz*12);
            
            % add interface nodes
            create_interface_nodes();
            B = B + T.L_interface*coef;
            
        end; % if compact
        
        build_dAdI();

        T.L = L+B;
        T.L = T.L/coef;
    end; % if (2==order)

    build_f0();
    
    function create_interface_nodes()

        l_iface_inds = (-3:3)*T.M;
        
        dgs_iface_t = ones(4*T.M,1) * T.L_z_interface_line;
        
        L_interface_loc1 = spdiags(dgs_iface_t, l_iface_inds, 7*T.M, 4*T.M).';
        
        L_interface_loc = L_interface_loc1(1+3*T.M:4*T.M,:);
        
        
        s1 = (T.interface_distance)*T.M;
        s2 = T.NM - 2*(s1+4*T.M);
        T.L_interface = blkdiag(...
            sparse(4*T.M,T.M), ...
            L_interface_loc, ...
            sparse((T.N-10)*T.M,(T.N-16)*T.M), ...
            L_interface_loc, ...
            sparse(4*T.M,T.M) );
        
    end; %create_interface_nodes()
end % function build_L()

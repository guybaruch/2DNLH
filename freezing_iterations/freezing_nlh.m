% Copyright (C) 2013, Guy Baruch.
% function freezing_nlh
%
% Solves tha NLH using the frosen Kerr fixed point iterations.
% at each iteration, the variable-coefficient Helmholtz is solved directly.
%    (real) K0, X_rp, Z_rp, xmax, ymax, p=N/Nc,
%    (int) max_I, max_J, flags=FR4S

function [xs zs E] = freezing_nlh(k0, Z_max, N, X_max, M, epsilon, SBC, order, compact, max_I )

    global T;
    
    addpath('../shared', '../shared/longitudinal', '../shared/transverse');
    
    T.radial = 0;
    T.epsilon = epsilon;
    T.sigma = 1;
       T.d_tau = 0.5;
    T.interface_E_inc_fix = 0;
    T.echo_args = 1;
    T.max_I=max_I;
    T.plot = 1;
    
    x_asymmetric = 0
    half_int_xs = 1

    init_grid(k0, Z_max, N, X_max, M, x_asymmetric, half_int_xs, ...
              SBC, order, compact)
    
    %generate_E_inc(k0, 'load');
    generate_E_inc(k0, 'sech');

    build_L();
    
    T.plotit=1;
    T.trace_jumps = 5;
    T.N_plots=max_I/T.trace_jumps;
    T.initial_guess = 'zero';
    
    prep_test();
    
    get_initial_E();

    T.A = T.L;

    %    T.res = (T.L * T.V_cont + T.f)*T.h_z;
    %    T.res_E = reshape(T.res,T.M,T.N);
    %    res_max = max(abs(T.res))

    
    for iter = 1:max_I
    
        T.iter = iter;

        T.V = reshape(T.E, T.NM, 1);
        update_A();

        %figure(11); plot(T.Zs,(abs( T.A*T.u -T.f0)))
        
        T.E_prev = T.E;

        %V_ext = T.A\T.f0;
        [l,u,p,q,r] = lu(T.A);
        V_new = q*( u\(l\( p*(r\T.f0) )) );
        
        T.dV = (V_new - T.V)*T.d_tau;
        
        T.V = T.V + T.dV;

        T.E = reshape(T.V, T.M, T.N);
        
        T.res_V = T.A*T.V-T.f0;
        
        if (1==T.trace_jumps || 1==rem(T.iter,T.trace_jumps) )
            T.res = reshape(T.res_V, T.M, T.N);
            T.dE = reshape(T.dV, T.M, T.N);
            T.iter
            do_trace();
    do_plots();
        end;
        
    end; % for T.iter
    
    
    if (max_I) 
        %        T.df(T.t_i,4:5) 
    end;

    if (T.trace_jumps==0)
        T.res = reshape(T.res_V, T.M, T.N);
        T.dE = reshape(T.dV, T.M, T.N);
    end;

    do_plots();

    xs = T.xs; zs = T.zs; E = T.E;

end % freezing_nlh

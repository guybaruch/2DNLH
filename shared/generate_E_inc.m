function generate_E_inc(k0, inc_shape)
% GENERATE_E_INC Generates the input fields to the NLH BVP. These fields then
% change the RHS - see BUILD_F0
    global T;

    T.k0 = k0;
    
    T.Einc0 = zeros(T.M,1);
    T.Einc1 = zeros(T.M,1);

    T.use_Snell = 1;
    
    switch (inc_shape)
        case 'load'
            generate_E_inc_load();
        case 'gaussian'
            generate_E_inc_gaussian();
        case 'gaussian_shift_angle'
            generate_E_inc_gaussian_shifted_angle();
        case 'sech'
            generate_E_inc_sech();
        case 'sech_shift'
            generate_E_inc_sech_shifted();
        case 'sech_shift_angle'
            generate_E_inc_sech_shifted_angle();
        case 'counter_propogating_shifted_sechs'
            generate_E_inc_counter_propogating_shifted_sechs();
        case 'counter_propogating_angle_sechs'
            generate_E_inc_counter_propogating_angle_sechs();
        otherwise
            error('generate_E_inc kind invalid');
    end;

    if (T.interface_E_inc_fix)
        disp('rough fix for interface');
        if ( strcmp(inc_shape,'sech_shift_angle') || ...
            strcmp(inc_shape,'gaussian_shift_angle'))
            disp('fix affected by angle')
            n1 = sqrt(1+(abs2(T.Einc0).^T.sigma)*T.epsilon);
            % Snell law per x
            angs = asin( n1*sin(T.IC_angle));
            ct = cos(angs);
            nt = sqrt((abs2(T.Einc0).^T.sigma)*T.epsilon+ct.^2);
            coef = (nt+ct)./(2*ct);
            T.Einc0 = T.Einc0.*coef;
        else
            coef = ( sqrt(abs(T.Einc0).^(2*T.sigma)*T.epsilon+1) +1 )/2;
            T.Einc0 = T.Einc0.*coef;
        end;
    else
        disp('no fix for interface');
    end;
    
    
    
    
    function generate_E_inc_gaussian()
    
        disp('gaussian incoming beam')
        T.Einc0 = exp(-T.xs.^2).';
        
    end; % generate_E_inc_gaussian()

    function generate_E_inc_gaussian_shifted_angle()

        T.IC_angle = atan( (T.IC_shift) / (T.Z_max/2));

        % fix for snell law
        disp(['gaussian incoming beam shifted by ' ...
              num2str(T.IC_shift) ...
              ' at angle pi*' num2str(T.IC_angle/pi) ])
        xs_shift = T.xs'-T.IC_shift;
            L_df = T.k0; k0 = T.k0; theta = T.IC_angle
            T.Einc0 = exp(...
                xs_shift*(-j*k0+1/k0)*sin(theta) ...
                -(xs_shift.^2).*exp( ...
                    (xs_shift*sin(theta)/k0).^2 ...
                )*cos(theta)^2 ...
            );
    end; % generate_E_inc_sech_shifted_angle()

    function generate_E_inc_sech()
        
        disp('sech incoming beam')
        coef =  T.k0 * sqrt(T.epsilon) ;
        T.Einc0 = ones(T.M,1)./( cosh(T.xs'/sqrt(2)) * coef)  ;

    end; % generate_E_inc_sech()

    function generate_E_inc_sech_shifted()

        % E_z_max  =
        % exp(-i*k0*sin(Theta)*(Y-half_dist))*(epsilon*k0^2)^
        %    (-1/2)./cosh(cos(Theta)*(Y-half_dist)/sqrt(2));
 
        disp(['sech incoming beam shifted by ' num2str(T.IC_shift)])
        A0 =  1/(T.k0 * sqrt(T.epsilon)) ;
        xs_shift = T.xs'-T.IC_shift;
        T.Einc0 = A0*ones(T.M,1)./cosh(xs_shift/sqrt(2))  ;

    end; % generate_E_inc_sech_shifted()

    function generate_E_inc_sech_shifted_angle()

        % E_z_max  =
        % exp(-i*k0*sin(Theta)*(Y-half_dist))*(epsilon*k0^2)^
        %    (-1/2)./cosh(cos(Theta)*(Y-half_dist)/sqrt(2));

        A0 =  1/(T.k0 * sqrt(T.epsilon)) ;

            
        T.IC_angle = atan( (T.IC_shift) / (T.Z_max/2));

        % fix for snell law
        disp(['sech incoming beam shifted by ' ...
              num2str(T.IC_shift) ...
              ' at angle pi*' num2str(T.IC_angle/pi) ])
        xs_shift = T.xs'-T.IC_shift;
        if (0==T.use_Snell)
            c1 = (T.k0+1/(4*T.k0))*sin(T.IC_angle);
            c2 = sqrt(1/2)*cos(T.IC_angle);
            T.Einc0 = A0*exp(-i*c1*xs_shift)./ cosh(c2*xs_shift)  ;
        else    
            n1 = sqrt(1+(abs2(T.Einc0).^T.sigma)*T.epsilon);
            % Snell law per x
            phis = -i*cumsum(n1)*T.h_x*T.k0*sin(T.IC_angle);
            angs = asin( n1*sin(T.IC_angle));
            c2 = sqrt(1/2)*cos(T.IC_angle);
            T.Einc0 = A0*exp(phis)./ cosh(c2*xs_shift)  ;
        end;
    end; % generate_E_inc_sech_shifted_angle()

    function generate_E_inc_counter_propogating_shifted_sechs()
        disp('counter-propogating (anti-parallel) sechs ')
        generate_E_inc_sech_shifted();
        T.Einc1 = T.Einc0(T.M:-1:1);
        %T.Einc0 = zeros(T.M,1);
    end; % generate_E_inc_counter_propogating_shifted_sechs

    function generate_E_inc_counter_propogating_angle_sechs()
        disp('counter-propogating (at angle) sechs ')
        generate_E_inc_sech_shifted_angle();
        T.Einc1 = T.Einc0;
        %T.Einc0 = zeros(T.M,1);
    end; % generate_E_inc_counter_propogating_angle_sechs
    
    function generate_E_inc_load()
        disp(' loading Einc0 and Einc1 ');

        re0 = load('re_Einc0.dat','ascii');
        im0 = load('im_Einc0.dat','ascii');
        if (length(re0)~=length(im0)) error(); end;
        %if (length(re0)~=length(T.Einc0) error(); end;
        T.Einc0 = re0+j*im0;

        re1 = load('re_Einc1.dat','ascii');
        im1 = load('im_Einc1.dat','ascii');
        if (length(re1)~=length(im1)) error(); end;
        %if (length(re1)~=length(T.Einc1) error(); end;
        T.Einc1 = re1+j*im1;
    end; % generate_E_inc_load


end % generate_E_inc()

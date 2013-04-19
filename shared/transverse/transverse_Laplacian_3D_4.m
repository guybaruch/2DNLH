function transverse_Laplacian_3D_4()
% TRANSVERSE_LAPLACIAN_2D_4  Creates the fourth-order transverse Laplacian on a
% Cylindrical grid, for both Dirichlet BCs and local Sommerfeld BCs.
	
	global T;

	if (T.half_int_xs ==0 ) 
		error('unimplemented') 
	end;

	% spdiags( [ a b c d e] , -2:2,...)
	% c1 d2 e3
	% b1 c2 d3 e4
	% a1 b2 c3 d4 e5
	% -  a2 b3 c4 d5 e6
	
	e = ones(T.M, 1);
	rhos = ( (1:T.M)-0.5 )' ;
	l = rhos.^-1;

	%remeber that Dxx_stencil is divided by (12 eta_x^2)

	so4 = [ ... % fourth order five point stencils
		[1 -8 0 8 -1]/12 ; ... % D^4_r
		[-1 16 -30 16 -1]/12 ]; % D^4_rr

	so2 = [ ... % second order five point stencils
		[0 -1 0 1 0]/2 ; ...  % D^2_r
		[0 1 -2 1 0]; ... % D^2_rr
		[-1 2 0 -2 1]/2 ; ... % D^2_rrr
		[1 -4 6 -4 1] ]; % D^2_rrrr

		l = -l;

	if (T.compact==1)
		Dxx_stencil= ... %(note that trans -> -1 ) ?
			( e*so4(2,:)  + l*so4(1,:) )*12 ...
			- (e*so2(2,:) + l*so2(1,:) ) *T.eta_z^2 ...
	- ( e*so2(4,:) + 2*l*so2(3,:) - l.^2*so2(2,:) + l.^3*so2(1,:) ) * T.h_r^2 ;
		
		T.L_x_dgs_t = line2dgs_t(Dxx_stencil);

		Dxx_laplacian_stencil =  -( e*so2(2,:) + l*so2(1,:) )/12;
		T.L_x_Laplacian_dgs_t = line2dgs_t(Dxx_laplacian_stencil);

		interface_stencil = 12*( ...
			e*( [0 0 T.eta_x^2 0 0]+so4(2,:) ) + l*so4(1,:) ) ; 
		T.L_x_dgs_interface_t = line2dgs_t(...
			interface_stencil)*T.h_r/(22*T.eta_x);
	else
		Dxx_line = 12*(e*so4(2,:) + l*so4(1,:)); % 1/r D^4_r
		T.L_x_dgs_t = line2dgs_t(Dxx_line);
	end;
	T.L_x = spdiags(T.L_x_dgs_t, -2:2, T.M, T.M).';
	%full(T.L_x)
	T.L_x = T.L_x / (12 * T.eta_x^2);



    function dgs_t = line2dgs_t(dgs_t)
		%if (length(stencil)~=5) error(); end;

		%dgs_t = dgs_t_1; %ones(T.M,1)*stencil;
	
		% remainder of stencil at R=0
		RemStencilAxis = [dgs_t(1,5:-1:4) ; [0 dgs_t(2,5)]];

		% remainder of stencil at R=Rmax
		RemStencilRmax = [ [dgs_t(T.M-1,1) 0 ]; ...
			dgs_t(T.M,2:-1:1) ];

	
		% symmetric BC at axis 
		if (T.half_int_xs ==0 ) 
			error('no int_xs radial');
		end;

		% u_0=u_1, u_-1=u_2
		B0 = RemStencilAxis * [0 1; 1 0];

		dgs_t(1,2:3) = dgs_t(1,2:3) + B0(1,2:-1:1);
		dgs_t(2,3:4) = dgs_t(2,3:4) + B0(2,2:-1:1);
		
		if (0==T.trans_SBC)
			
			a = [0 -1 9 9 -1];  % u(n+1/2)=0
	
		else
	
			% Sommerfeld BC
			alpha = -T.k0 * besselh(1,1,T.k0*T.X_max) / besselh(0,1,T.k0*T.X_max) ;
			a = [0 1 -27 27 -1] - 1.5*alpha*T.h_x*[0 -1 9 9 -1] ;
			
		end; % SBC
		s = a(4)+4*a(5);
		v1 = [a(1)-a(5)   a(2)+4*a(5)   a(3)-6*a(5)] / s ;
		v2 = [4*a(1)+a(4) 4*a(2)-4*a(4) 4*a(3)+6*a(4)] / s ;
	
		M = T.M;
		B1 = - RemStencilRmax * [ v1 ; v2 ] ;
		%B1 = - [ ((M-0.5)/(M-1.5))*v1 ; ...
		%	     -(16*M/(M-0.5)) * v1 + ((M+0.5)/(M-0.5))*v2 ];
	
		% e d c b a
		% X e d c b8
	 	% X X e d c9
		
		% remove elements "outside matrix"
		dgs_t(1,4:5)=0;
		dgs_t(2,5)=0;
		dgs_t(T.M,1:2)=0;
		dgs_t(T.M-1,1)=0;
		
		
		%T.L_x(T.M-1:T.M,T.M-2:T.M) = [16 -30 16;-1 16 -30] + T.L_x_B1;
		dgs_t(T.M-1,2:4) = dgs_t(T.M-1,2:4) + B1(1,3:-1:1);
		dgs_t(T.M,3:5) = dgs_t(T.M,3:5) + B1(2,3:-1:1);
	end;

	
end % transverse_Laplacian_3D_4();

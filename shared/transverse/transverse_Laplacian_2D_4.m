function transverse_Laplacian_2D_4()
% TRANSVERSE_LAPLACIAN_2D_4  Creates the fourth-order transverse Laplacian on a
% Cartesian grid, for both Dirichlet BCs and local Sommerfeld BCs.
	
	global T;
	
	
	% spdiags( [ a b c d e] , -2:2,...)
	% c1 d2 e3
	% b1 c2 d3 e4
	% a1 b2 c3 d4 e5
	% -  a2 b3 c4 d5 e6
	
	e = ones(T.M, 1);

	if (T.compact==1)
		Dxx_line = ...
			[-1 16 -30 16 -1] ... % 12 D^4_xx
			- [1 -4 6 -4 1]*T.h_r^2 ... % h_z^2 D^2_xxxx
			- [0 1 -2 1 0] * T.eta_z^2; % eta_z^2 D^2_xx 

		Dxx_laplacian_line = - [0 1 -2 1 0]/12;
		T.L_x_Laplacian_dgs_t = line2dgs_t(Dxx_laplacian_line);

		%interface_line = [ 0 1 -2+T.eta_x^2 1 0 ] ;  % O(h^3) BC
		interface_line = [-1 16 -30+12*T.eta_x^2 16 -1]; % O(h^4) BC
		T.L_x_dgs_interface_t = line2dgs_t(interface_line)*T.h_r/(22*T.eta_x);

	else
		Dxx_line = [-1 16 -30 16 -1];
	end;
	T.L_x_dgs_t = line2dgs_t(Dxx_line);
	T.L_x = spdiags(T.L_x_dgs_t, -2:2, T.M, T.M).';
	%full(T.L_x);
	T.L_x = T.L_x / (12 * T.eta_x^2);



	
    function dgs_t = line2dgs_t(stencil)
		if (length(stencil)~=5) error(); end;

		dgs_t = ones(T.M,1)*stencil;
		% remove elements "outside matrix"
		dgs_t(1,4:5)=0;
		dgs_t(2,5)=0;
		dgs_t(T.M,1:2)=0;
		dgs_t(T.M-1,1)=0;
	
		% [u_M]
		
		if (0==T.trans_SBC)
			if (T.half_int_xs )
				% u_M=-u_M-1, u_M+1=-u_M-2
				C = [0 0 -1; 0 -1 0 ];
			else
				% u_M=0, u_M+1=-u_M-1
				C = [0 0 0; 0 0 -1 ];			
			end;
			C5 = 0; C6 = 0;
			C2 = C(1,2); 	C1 = C(1,3);
			C4 = C(2,2); 	C3 = C(2,3);
			
		else
			if (T.half_int_xs ==0 ) 
				error('unimplemented') 
			end;
	
			% Sommerfeld BC
			q= [0 0];
			%d1 = 8-sqrt(36+12*T.eta_x^2); q(1) = (d1+i*sqrt(4-d1^2))/2;
			d1 = 8-sqrt(36+12*T.eta_x^2); q(1) = (d1+sqrt(d1^2-4))/2;
			d2 = 8+sqrt(36+12*T.eta_x^2); q(2) = (d2-sqrt(d2^2-4))/2;
			
			%		q(1) = exp(j*T.eta_x);
			
			C5 = 0; C6 = 0;
			BC_type = 3;
			switch (BC_type)
				case 0
					% guybar: this is from FT04,
					s = q(2)+q(2)^2-q(1)*(1+q(2)^3);
					C1 = (q(2)+q(2)^2-q(1)^2*(1+q(2)^3))/s;
					C2 = -(q(2)+q(2)^2)*(q(1)-q(1)^2)/s;
					C3 = (1+q(2)^3)*(1-q(1)^3)/s;
					C4 = -(q(1)*(1+q(2)^3)-q(1)^3*(q(2)+q(2)^2))/s;
				case 1
					% also from FT04 : older approach
					C1 = q(1)+q(2) ;
					C2 = -q(1)*q(2);
					C3 = (q(1)+q(2))^2-q(1)*q(2);
					C4 =  -(q(1)+q(2))*q(1)*q(2);
				case 2
					% as in BFT:06, third order
					c = [1 -27 27 -1] - 1.5*T.eta_x*j* [-1 9 9 -1]; % -2->1
	
					s = -1/( c(3)+3*c(4) )	;
	
					C1 = ( c(2)-3*c(4) ) *s;
					C2 = ( c(1)+c(4) ) *s;
					C3 = ( 3*(c(2)+c(3)) ) *s;
					C4 = ( 3*c(1)-c(3) ) *s;
				case 3
					% as in BFT:06, fourth order
					c = [1 -27 27 -1] - 1.5*T.eta_x*j* [-1 9 9 -1]; % -2->1
	
					s = -1/( c(3)+4*c(4) )	;
	
					C5 = ( -c(4) ) *s;
					C6 = ( c(3) ) *s;
					C2 = ( c(1)+4*c(4) ) *s;
					C1 = ( c(2)-6*c(4) ) *s;
					C4 = ( 4*(c(1)-c(3)) ) *s;
					C3 = ( 4*c(2)+6*c(3) ) *s;
				otherwise
					error('unimplemented');
			end; % switch
			
			C = [C5 C2 C1; C6 C4 C3];
			
		end; % SBC
	
		%B1 = [-1 0; 16 -1]*C;
		% remainder of stencil at Xmax
		RemStencilXmax = [[stencil(5) 0]; stencil(4:5)];
		B1 = RemStencilXmax*C;
		
		% e d c b a
		% X e d c b8
	 	% X X e d c9
		
		%T.L_x(T.M-1:T.M,T.M-2:T.M) = [16 -30 16;-1 16 -30] + T.L_x_B1;
	%	dgs_t(T.M-1:T.M,3:5) = dgs_t(T.M-1:T.M,3:5) + B1(:,1:3);
		dgs_t(T.M-1,2:4) = dgs_t(T.M-1,2:4) + B1(1,3:-1:1);
		dgs_t(T.M,3:5) = dgs_t(T.M,3:5) + B1(2,3:-1:1);
		
		% remainder of stencil at Xmin
		% RemStencilXmin = [-1 16; 0 -1];
		RemStencilXmin = [stencil(1:2) ; [0 stencil(1)]];
	
		% BC at lower x
		if (T.x_asymmetric)
			% same as X_min, just reordered
			C = [C3 C4 C6; C1 C2 C5];
		else
			% symmetric BC at axis 
			if (T.half_int_xs ) 
				% u_0=u_1, u_-1=u_2
				C = [0 1 0; 1 0 0];
			else
				% u_0=u_2, u_-1=u_3
				C = [0 0 1; 0 1 0];
			end;
		end;
		B0 = RemStencilXmin*C;
	
		%T.L_x(1:2,1:3) = T.L_x(1:2,1:3)+T.L_x_B0; 
		dgs_t(1,1:3) = dgs_t(1,1:3) + B0(1,3:-1:1);
		dgs_t(2,2:4) = dgs_t(2,2:4) + B0(2,3:-1:1);
	
	end; % line2dgs_t
	
end % transverse_Laplacian_2D_4();

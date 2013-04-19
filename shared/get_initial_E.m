function get_initial_E()
% constructs the initial guess.
% T.initial_guess = 'zero', 'keep', 'load'.
	
	global T;
	
	%T.E0 = zeros(T.M,T.N);
	
	switch (T.initial_guess)
		case 'zero'
			disp('E0 = 0');
			T.E = zeros(T.M,T.N);
		case 'keep'
			disp('keeping E');
		case 'load'
			disp('');
			rE1 = load('./old/re_E.dat','-ascii');
			iE1 = load('./old/im_E.dat','-ascii');
			E1 = rE1+i * iE1;

			xs1 = load('./old/xs.dat','-ascii');
			zs1 = load('./old/zs.dat','-ascii');
			M1 = length(xs1);
			N1 = length(zs1);

			display(sprintf('loading E0 using interp from a field with (N,M)=(%d,%d)',N1,M1));
			
			%imethod = '*linear';
			imethod = '*spline';
			%imethod = '*cubic';

						
			dx1 = xs1(2)-xs1(1);
			xs1 = [xs1 ; xs1(M1) + dx1*2];
			E1_xmax = 3*E1(M1,:) -2*E1(M1-1,:);
			E1 = [E1 ; E1_xmax];

			if (xs1(1)>0)
				xs1 = [0; xs1]; 
				E1_x0 = (9*E1(1,:) -E1(2,:) )/8;
				E1 = [E1_x0 ; E1 ];
			end;
			
			
			dz1 = zs1(2)-zs1(1);
			zs1 = [zs1(1)-2*dz1 ; zs1 ; zs1(N1) + 2*dz1 ];
			E1_z0 = 3*E1(:,1) -2*E1(:,2);
			E1_zmax = 3*E1(:,N1) -2*E1(:,N1-1);
			E1 = [E1_zmax  E1  E1_z0];
			
			%[size(E1) size(xs1) size(zs1)]
			
			
			T.E0 = interp2(zs1,xs1, E1, T.zs, T.xs',imethod,0);
			T.E = T.E0;
		case 'NLS'
			disp('using coupled NLS for initial Guess');
			get_initial_E_NLSises();
		otherwise
			error('unimplemented');
	end;

end %get_initial_E;
	

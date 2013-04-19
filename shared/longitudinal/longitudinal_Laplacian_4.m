function longitudinal_Laplacian_4()
% LONGITUDINAL_LAPLACIAN_4  Creates the fourth-order 1D Laplacian L_z
% using the standard [-1 16 -30 16 -1] centered-difference stencil.

	global T;
	
	T.L_z = [-1 16 -30 16 -1] / (12*T.eta_z^2);

	e = ones(T.M,1);


	% Sommerfeld BC
	q= [0 0];
	%[ T.eta_pars.^2 8-sqrt(36+12*T.eta_pars.^2) ]

	d1 = 8-sqrt(36+12*T.eta_pars.^2);    q1 = (d1+i*sqrt(4-d1.^2))/2;
	% for very negative k_parallel, i.e. LARGE negative eta_pars.^2, 
	% we have imag(d1)<0  and so the "wrong" root is chosen
	% im(sqrt(d1^2-4))<0, while im(i*sqrt(4-d1^2))>0
 	
	
	d2 = 8+sqrt(36+12*T.eta_pars.^2);    q2 = (d2-sqrt(d2.^2-4))/2;
	
	if (0)
		C1 = q1+q2;
		C2 = -q1.*q2;
		C3 = (q1+q2).^2 - q1.*q2;
		C4 = -(q1+q2).*q1.*q2;
	else
		% guybar: this is from FT04,
		C1 = (q2+q2.^2-q1.^2.*(1+q2.^3))./(q2+q2.^2-q1.*(1+q2.^3));
		C2 = -(q2+q2.^2).*(q1-q1.^2)./(q2+q2.^2-q1.*(1+q2.^3));
		C3 = (1+q2.^3).*(1-q1.^3)./(q2+q2.^2-q1.*(1+q2.^3));
		C4 = -(q1.*(1+q2.^3)-q1.^3.*(q2+q2.^2))./(q2+q2.^2-q1.*(1+q2.^3));
	end;

	T.G1 = diag(C1);
	T.G2 = diag(C2);
	T.G3 = diag(C3);
	T.G4 = diag(C4);

	%sq1 = sqrt(q1);
	%sq1 = ones(T.M,1)./sqrt(q1)
	%v = conj([ sq1 sq1.*q1 sq1.*(q1.^2) sq1.*(q1.^3) ])
	
	sq1 = sqrt(q1);
	v = [ sq1 sq1.*q1 sq1.*(q1.^2) sq1.*(q1.^3) ];
	v = ones(T.M,4)./v;

	
	phi_0 = C4.*v(:,1) + C3.*v(:,2) - v(:,4);
	phi_1 = C2.*v(:,1) + C1.*v(:,2) - v(:,3);

	T.F0 = diag(-phi_0);
	T.F1 = diag(-phi_1);

	z = zeros(T.M);
	e = eye(T.M);
	
	T.Bz0 = [-e 16*e; z -e] * ...
			[ U2E(T.G3) U2E(T.G4); U2E(T.G1) U2E(T.G2) ];
	T.Bf0 = [-e 16*e; z -e] * ...
			[ U2E(T.F0) z; z U2E(T.F1) ] ;

	T.Bz1 = [-e z; 16*e -e] * ...
			[ U2E(T.G2) U2E(T.G1); U2E(T.G4) U2E(T.G3) ];
	T.Bf1 = [-e z; 16*e -e] * ...
			[ U2E(T.F1) z; z U2E(T.F0) ];

end

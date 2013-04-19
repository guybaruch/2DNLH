function build_f0()
% BUILD_F0  Given the incoming fields Einc0 and Einc1, builds the RHS f0
% in A(E)E-f0=0.
	global T;
	
	if (2==T.order || T.compact)
		T.f0(1:T.M) = -T.Bf * T.Einc0 / T.eta_z^2;
		T.f0(T.NM-T.M+1:T.NM) = -T.Bf * T.Einc1 / T.eta_z^2;
	else
		T.f0(1:2*T.M) = -T.Bf0 * [T.Einc0; T.Einc0] / (12*T.eta_z^2);
		T.f0(T.NM-2*T.M+1:T.NM) = -T.Bf1 * [T.Einc1; T.Einc1] / (12*T.eta_z^2);
	end;

end % build_f0()

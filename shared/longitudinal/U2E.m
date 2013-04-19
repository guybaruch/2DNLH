function E = U2E(U)
% U2E  Fourier tranforms the field in the transverse direction - see the
% definition of the longitudinal BCs in the 2009 JCP.

	global T;
	% Psi F Psi^-1
%	E = T.Psi_U\(T.Psi_L\(T.Psi_P * U));
	E = T.Psi*U*T.Psi_inv;
end


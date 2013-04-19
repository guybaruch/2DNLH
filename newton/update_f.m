function update_f()
% UPDATE_F  updates the RHS of the NLH functional Av=f
    global T;
	F = T.A*T.V-T.f0;
	T.f_2 = c2r(F);
end % update_f

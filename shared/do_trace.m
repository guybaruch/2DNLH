function do_trace() 
% DO_TRACE  Performs the traces in the course of an iterative solver run.
	global T;

    if (0==T.trace_jumps) return; end;
    
	T.t_i = T.t_i+1;

	T.df(T.t_i,1)=max(abs(T.dV));
	T.daE=abs(T.E)-abs(T.E_prev);
	T.df(T.t_i,2)=max(max(abs(T.daE)));
	
	T.res_V=(T.A*T.V-T.f0)*T.h_z;
	T.df(T.t_i,3)=sqrt(max(abs2(T.res_V)));

	T.df(T.t_i,4)=sqrt(max(abs2(T.V)));

	if (T.plot==2)
		do_plots();
	end;
	
	T.df(T.t_i,:)
	
end % do_trace()

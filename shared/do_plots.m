% function do_plots()
% 
%	plots results of 1D run, whether Newton or fixed-point. 
%	Assumes T is a global variable

function do_plots()
	
	plot_Poynting = 1;
	plot_trace = 0;
	plot_last_diff = 0;
	plot_on_axis_amp = 0;
	plot_residual = 1;
	
	global T;

		
	switch (T.plot)
		case 0
			return
		case 1
			if (T.iter < T.max_I)
				return
			end;
	end;
	
	% choose the grid-scale
	switch T.normalize_grid
		case 0
			gs = 1; 
		case 1
			gs = T.Lambda0;
		case 2
			gs = T.k0*2;
		otherwise
			error('');
	end;

	zs1 = T.zs/T.gs;
	xs1 = T.xs/T.gs;
	
	figure(1);  
	surfc(zs1,xs1,abs(T.E),'EdgeColor', 'none')
	title('|E|');


	if (plot_last_diff)
		figure(7);  
		surf(zs1,xs1,abs(T.dE),'EdgeColor', 'none')
		title('last diff');

		%t = [Zs'  T.daE]; save 'run/last_diff.dat' -ascii t;
	end;
	

	if (plot_residual)
		figure(6);
		ar = abs(T.res); 
		surf(zs1,xs1,ar,'EdgeColor', 'none');
		title('residual');
	end;
	
	%t = [Zs' ar]; save 'run/power_res.dat' -ascii t;

	if plot_on_axis_amp
		figure(5);
		plot(zs1, abs(T.E(1,:)),'r');
		hold on; plot(zs1, abs(T.E(2,:)),'b');
		title('On axis amplitude');
		hold off;
	end;
	
	N_trace = length(T.df(:,1));
	%T.Nz = length(T.u);
	
	if (N_trace>2 && plot_trace)
		figure(4);
		t_is = 1:N_trace;
		t_is = t_is * T.trace_jumps;
		semilogy(t_is, T.df(:,1),'r');
		hold on;
		semilogy(t_is, T.df(:,2),'r--');
		semilogy(t_is, T.df(:,3),'b');
		%		semilogy(t_is, T.df(:,4),'k');
		%		semilogy(t_is, T.df(:,5),'k--');
		legend('dE','d|E|','residual')	
		title('convergence');
		hold off;
	end;

	if (plot_Poynting)
		[T.Ez T.Ex T.Sz T.Sx T.curlS T.divS] = Poyntings(T.E, T.h_z, T.h_x);
		
		figure(2); surfc(zs1, xs1, T.Sz, 'EdgeColor','none'); 
		title('Sz'); %ylim([0 2]); 
		figure(3); contourf(zs1, xs1, T.Sx); %ylim([0 2]);
		title('Sx'); colorbar; 
	end;
	
		
end % do_plots()

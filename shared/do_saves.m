function do_saves()
% DO_SAVES  Save the current grid and fields as ascii files, separated to real
% and imaginary components.
    global T;
    
    if (T.save)
    
        disp('saving ...');
        
        t = T.xs'; save_it('run/xs.dat', 1);
        t = T.zs'; save_it('run/zs.dat', 1);


        t = T.xs'/T.gs; save_it('run/xs-Ldf.dat', 1);
        t = T.zs'/T.gs; save_it('run/zs-Ldf.dat', 1);
        
        t = real(T.E); save_it( 'run/re_E.dat', T.save_double);
        t = imag(T.E); save_it( 'run/im_E.dat', T.save_double);
        t = abs(T.E); save_it( 'run/abs_E.dat', T.save_double);
        t = abs(T.dE); save_it( 'run/abs_dE.dat', T.save_double);
    
        if (T.half_int_xs)
            axis_amp = abs(9*T.E(1,:)-T.E(2,:))/8;
            axis_Sz = (9*T.Sz(1,:)-T.Sz(2,:))/8;
        else
            axis_amp = abs(T.E(1,:));
            axis_Sz = T.Sz(1,:);
        end;
        t = [T.zs/T.gs; axis_amp]'; save_it( 'run/axis_amp.dat', 1 );
        t = [T.zs/T.gs; axis_amp.^2 ; axis_Sz/T.k0]'; save_it( 'run/axis_amp2_Sz.dat', 1 );
        if (T.trace_jumps)
            t = T.df; save_it( 'run/dfs', 0);
        end;
    end;
    
    function save_it(name, dbl)
        %disp([ 'saving ' name]);
        if (dbl)
            save(name,'t','-ascii','-double');
        else
            save(name,'t','-ascii');
        end;
    end;
end

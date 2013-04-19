
lambda0 = 2*pi/k0

if radial
	r0 = 1
	L_df = k0
else
	r0 = sqrt(2)
	L_df = 2*k0
end;

% close all; clear all; pack; radial=1; sym_field=1; epsilon=0.15; k0=8; Z_max=9; X_max = 3; curl_div =1; fft_center = 1; fwd_bkwd=1;

rE1 = load('./old/re_E.dat','-ascii');
iE1 = load('./old/im_E.dat','-ascii');

E1 = rE1+i * iE1;

xs = load('./old/xs.dat','-ascii');
zs = load('./old/zs.dat','-ascii');
M = length(xs)
N = length(zs)

cut_anodes = 1

if (cut_anodes)
	if compact
		E1 = E1(:,5:N-4);
		zs = zs(5:N-4);
		N=N-8
	else
		E1 = E1(:,3:N-2);
		zs = zs(3:N-2);
		N=N-4
	end
end

if (x_asymmetric==0)
	E = zeros(M*2+1,N);
	E(1:M,:) = E1(M:-1:1,:);
	E(M+1,:) = (9*E1(1,:)-E1(2,:))/8;
	E(M+2:2*M+1,:) = E1;
	xs =  [ -xs(M:-1:1) ; 0 ; xs];
%	M = 2*M;
else
	E = E1;
end;



h_x = xs(2)-xs(1)
h_z = zs(2)-zs(1)

[Ez Ex Sz Sx curlS divS] = Poyntings(E, h_z, h_x);
%
%if (x_asymmetric==0)
%	figure(1); 
%	subplot(1,2,1);
%	plot(zs, Sz(M+1,:)); hold on; plot(zs, Sz(M+2,:),'r'); plot(zs, Sz(M+3,:),'g'); hold off; title('Sz on axis');
%	subplot(1,2,2); plot(zs, Sx(M+1,:)); hold on; plot(zs, Sx(M+2,:),'r'); plot(zs, Sx(M+3,:),'g'); hold off; title('Sx on axis');
%end;
%
%
%
%figure(3);  surfc(zs, xs, abs(E),'EdgeColor','none'); 
%title('amplitude |E|','FontSize',18); set(gca,'FontSize',16); xlim([0,Z_max])
%%xlabel('Z/L_{DF}','FontSize',20); ylabel('X/\lambda_{0}','FontSize',20,'Rotation',0); set(gca,'FontSize',16);
%figure(4); surfc(zs, xs, Sz,'EdgeColor','none'); 
%title('forward Poynting S_Z','FontSize',18); set(gca,'FontSize',16); xlim([0,Z_max])
%figure(5); surfc(zs, xs, Sx,'EdgeColor','none');  title('transverse Poynting Sx')
%
%%if exists(curl_div)
%	%figure(6); surfc(zs, xs, curlS,'EdgeColor','none'); title('curl S');
%	%figure(7); surfc(zs, xs, divS,'EdgeColor','none'); title('div S');
%%end;
%
%if (radial)
%	t = ones(M,N).*(abs(xs(1:M))*ones(1,N));
%	Nc = 7.4492/ (epsilon*k0^2)  % 1.8623
%	Nc = 1.8623/(epsilon*k0^2) % 1.8623
%else
%	t = ones(M,N);
%	Nc = 1
%end;
%
%Sz1 = Sz(1:M,:).*t*h_x/k0;
%
%tot_pow = sum(Sz1,1);
%core_pow = sum(Sz1(0.8*M:M,:),1);
%figure(8); plot(zs, tot_pow/Nc,'b'); xlim([0 1])
%hold on; plot(zs, core_pow/Nc,'r'); 
%title('total and core \int Sz dx'); hold off;
%
%%a = abs(E(M,N-30:N)); b = zs(N-30:N)-zs(N); c = [a' b];
%%figure(8); hold on;  plot(b, a)
%%save 'oaE2_suf_920.dat' -ascii c
%
%if (0)
%	ks = (-N/2:(N-2)/2)*2*pi/Z_max;
%	v1 = fftshift(fft((E)));
%	t = [ks'/k0 abs(v1)'];
%	save 'axis_E_fft.dat' -ascii t
%	figure(10); semilogy(ks/k0, abs(v1));
%end;
%
%if (0)
%
%	E_fwd = (E+Ez/(k0*j))/2;
%	E_bkwd = (E-Ez/(k0*j))/2;
%	%figure(13); surfc(zs, xs, abs(E_fwd),'EdgeColor','none'); title('forward moving |E_{fwd}|')
%	%figure(14); surfc(zs, xs, abs(E_bkwd),'EdgeColor','none'); title('backward moving |E_{bkwd}|')
%	
%	wf_step = N/20
%	figure(15); waterfall(xs,zs(1:wf_step:N), abs(E_fwd(:,1:wf_step:N))');
%	title('forward moving |E_{fwd}|')
%	figure(16); waterfall(xs,zs(1:wf_step:N), abs(E_bkwd(:,1:wf_step:N))');
%	title('backward moving |E_{bkwd}|')
%
%end
%
%if (1)
%	figure(20); plot(xs,Sz(:,1),'b',xs,Sz(:,N),'r');
%	%xlim([-10,10]); ylim([-1.5,1.5]);
%	set(gca,'FontSize',16); title('Forward Poynting S_Z at interfaces','FontSize',18)
%end
%

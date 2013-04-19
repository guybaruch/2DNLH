function [Ez Ex Sz Sx curlS divS] = Poyntings(E, h_z, h_x)
	[Ez Ex] = gradient(E, h_z, h_x);

	Sx = imag(Ex.*conj(E));
	Sz = imag(Ez.*conj(E));
	
	[Sxz Sxx] = gradient(Sx, h_z, h_x);
	[Szz Szx] = gradient(Sz, h_z, h_x);
	
	curlS = Sxz - Szx;
	
	divS = Sxx + Szz;
	
end
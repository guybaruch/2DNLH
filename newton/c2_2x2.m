function m2x2 = c2_2x2(c)
% C2_2X2  Multiplication by complex number represented in 2x2 matrix form
%
%	I_2 = [1 0 ; 0 1 ];
%	i_sigma2 = [ 0 1 ; -1 0 ];
%	m2x2 = (real(c) * I_2) - (imag(c) * i_sigma2);
	
%	m2x2 = [ real(c) -imag(c); imag(c) real(c)];
	re = real(c); im = imag(c);
	m2x2 = [ re -im; im re];
end % c2_2x2

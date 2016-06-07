function te = expand(t,f)

% Expand spatially an image t in a factor f 
% in X and in Y.
% t may be complex.
% It fills in with zeros in the Fourier domain.
%	te = expand(t, f)
% See also: shrink.m
% JPM, May 95, Instituto de Optica, CSIC, Madrid.

[my mx]=size(t);
my=f*my;
mx=f*mx;
Te=zeros(my,mx);
T=f^2*fftshift(fft2(t));
y1=my/2+2-my/(2*f);
y2=my/2+my/(2*f);
x1=mx/2+2-mx/(2*f);
x2=mx/2+mx/(2*f);
Te(y1:y2,x1:x2)=T(2:my/f,2:mx/f);
Te(y1-1,x1:x2)=T(1,2:mx/f)/2;
Te(y2+1,x1:x2)=((T(1,mx/f:-1:2)/2)').';
Te(y1:y2,x1-1)=T(2:my/f,1)/2;
Te(y1:y2,x2+1)=((T(my/f:-1:2,1)/2)').';
esq=T(1,1)/4;
Te(y1-1,x1-1)=esq;
Te(y1-1,x2+1)=esq;
Te(y2+1,x1-1)=esq;
Te(y2+1,x2+1)=esq;
Te=fftshift(Te);
te=ifft2(Te);
if all(imag(t)==0),
	te = real(te);
end

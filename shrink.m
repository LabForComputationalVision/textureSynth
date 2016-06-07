function ts=shrink(t,f)

% It shrinks an image in a factor f
% in each dimension. 
%	ts = shrink(t,f)
% ts may also be complex.
% See also: expand.m, blurDn.m
% JPM, April 95, Instituto de Optica, CSIC, Madrid.

[my,mx]=size(t);
T=fftshift(fft2(t))/f^2;
Ts=zeros(my/f,mx/f);
y1=my/2+2-my/(2*f);
y2=my/2+my/(2*f);
x1=mx/2+2-mx/(2*f);
x2=mx/2+mx/(2*f);
Ts(2:my/f,2:mx/f)=T(y1:y2,x1:x2);
Ts(1,2:mx/f)=(T(y1-1,x1:x2)+T(y2+1,x1:x2))/2;
Ts(2:my/f,1)=(T(y1:y2,x1-1)+T(y1:y2,x2+1))/2;
Ts(1,1)=(T(y1-1,x1-1)+T(y1-1,x2+1)+T(y2+1,x1-1)+T(y2+1,x2+1))/4;
Ts=fftshift(Ts);
ts=ifft2(Ts);
if all(imag(t)==0),
	ts = real(ts);
end

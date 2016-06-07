function [Y,snrV,Chf]=modacor22(X,Cy,p);

% It imposes the desired autocorrelation in the given (central) samples (Cy) to
% an image X, convolving it with an even filter of size(Cy), in such a way
% that the image containts change as less as possible, in a LSE sense.
%	[Y,snr,Chf]=modacor22(X,Cy,p);
%	Chf:	Fourier transform of the filter that forces the autocorrelation
%	p [OPTIONAL]:	mixing proportion between Cx and Cy
%			it imposes (1-p)*Cx + p*Cy,
%			being Cx the actual autocorrelation.
%			DEFAULT: p = 1; 

% JPM, 10/97, working with EPS, NYU

Warn = 0;  % Set to 1 if you want to see warning messages
if (exist('p') ~= 1)
  p = 1;
end

% Compute the autocorrelation function of the original image

[Ny,Nx]=size(X);
Nc=size(Cy,1); 	% Normally Nc<<Nx, only the low indices of the autocorrelation	
if (2*Nc-1 > Nx) & Warn
  warning('Autocorrelation neighborhood too large for image: reducing');
  Nc = 2*floor(Nx/4)-1;
  first = (size(Cy,1)-Nc)/2;
  Cy = Cy(first+1:first+Nc, first+1:first+Nc);
end

Xf=fft2(X);
Xf2=abs(Xf).^2;
Cx=fftshift(real(ifft2(Xf2)))/(2-isreal(X));
Cy=Cy*prod(size(X));	% Unnormalize the previously normalized correlation

cy=Ny/2+1;
cx=Nx/2+1;
Lc=(Nc-1)/2;
Cy0 = Cy;
Cy = p*Cy + (1-p)*Cx(cy-Lc:cy+Lc,cx-Lc:cx+Lc);

% Compare the actual correlation with the desired one
%imStats(Cx(cy-Lc:cy+Lc,cx-Lc:cx+Lc),Cy)
snrV=10*log10(sum(sum(Cy0.^2))/sum(sum((Cy0-Cx(cy-Lc:cy+Lc,cx-Lc:cx+Lc)).^2)));

% Take just the part that has influence on the samples of Cy (Cy=conv(Cx,Ch))
Cx=Cx(cy-2*Lc:cy+2*Lc,cx-2*Lc:cx+2*Lc);

% Build the matrix that performs the convolution Cy1=Tcx*Ch1

Ncx=4*Lc+1;
M=(Nc^2+1)/2;
Tcx=zeros(M);

for i=Lc+1:2*Lc,
	for j=Lc+1:3*Lc+1,
		nm=(i-Lc-1)*(2*Lc+1)+j-Lc;
		ccx=Cx(i-Lc:i+Lc,j-Lc:j+Lc);
		ccxi=ccx(2*Lc+1:-1:1,2*Lc+1:-1:1);
		ccx=ccx+ccxi;
		ccx(Lc+1,Lc+1)=ccx(Lc+1,Lc+1)/2;
		ccx=vectify(ccx');
		Tcx(nm,:)=ccx(1:M)';
	end
end
i=2*Lc+1;
for j=Lc+1:2*Lc+1,
	nm=(i-Lc-1)*(2*Lc+1)+j-Lc;
	ccx=Cx(i-Lc:i+Lc,j-Lc:j+Lc);
	ccxi=ccx(2*Lc+1:-1:1,2*Lc+1:-1:1);
	ccx=ccx+ccxi;
	ccx(Lc+1,Lc+1)=ccx(Lc+1,Lc+1)/2;
	ccx=vectify(ccx');
	Tcx(nm,:)=ccx(1:M)';
end

% Rearrange Cy indices and solve the equation

Cy1=vectify(Cy');
Cy1=Cy1(1:M);

Ch1=inv(Tcx)*Cy1;

% Rearrange Ch1

Ch1=[Ch1;Ch1(length(Cy1)-1:-1:1)];
Ch=reshape(Ch1,Nc,Nc)';

% Compute H from Ch (H is zero-phase) through the DFT

%s=2^(ceil(log(Nc)/log(2))+1);
%H=sqrt(abs(fft2(Ch,s,s)));
%h=fftshift(real(ifft2(H)));
%h=h(s/2+1-Lc:s/2+1+Lc,s/2+1-Lc:s/2+1+Lc);
%%plot(Ch);drawnow
%h=recphase(Ch);

% Compute Y as conv(X,H) in the Fourier domain

%%Y=real(ifft2(Xf.*H));
%Y=real(ifft2(Xf.*sqrt(abs(fft2(Ch,Ny,Nx)))));
aux=zeros(Ny,Nx);
aux(cy-Lc:cy+Lc,cx-Lc:cx+Lc)=Ch;
Ch=fftshift(aux);
Chf=real(fft2(Ch));
%Chf=fft2(Ch,Ny,Nx);
%figure(7);plot(Chf);drawnow;
Yf=Xf.*sqrt(abs(Chf));
Y=ifft2(Yf);
%Y=cconv2(X,h);

% Checks the fidelity of the imposition

%Cy2=fftshift(real(ifft2(Xf2.*abs(Chf))))/(2-isreal(X));
%Cy2=Cy2(cy-Lc:cy+Lc,cx-Lc:cx+Lc);
%imStats(Cy,Cy2)
%imStats(X,Y)

%Yf2=abs(Yf).^2;
%Cy3=fftshift(real(ifft2(Yf2)))/2;

%Cy3=Cy3(cy-Lc:cy+Lc,cx-Lc:cx+Lc);
%snr(Cy,Cy-Cy3)
%imStats(Cy,Cy3)


		

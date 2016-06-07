function [chm, snrk] = modskew(ch,sk,p);

% Adjust the sample skewness of a vector/matrix, using gradient projection,
% without affecting its sample mean and variance.
%
% This operation is not an orthogonal projection, but the projection angle is
% near pi/2 when sk is close to the original skewness, which is a realistic
% assumption when doing iterative projections in a pyramid, for example
% (small corrections to the channels' statistics).
%
% 	[xm, snrk] = modskew(x,sk,p);
%		sk: new skweness
%       	p [OPTIONAL]:   mixing proportion between sk0 and sk 
%               	        it imposes (1-p)*sk0 + p*sk,
%                      		being sk0 the current skewness.
%                       	DEFAULT: p = 1;

%
% JPM. 2/98, IODV, CSIC
% 4/00, CNS, NYU

Warn = 0;  % Set to 1 if you want to see warning messages
if ~exist('p'), 
  p = 1;
end

N=prod(size(ch));	% number of samples
me=mean2(ch);
ch=ch-me;

for n=2:6,
	m(n)=mean2(ch.^n);
end

sd=sqrt(m(2));	% standard deviation
s=m(3)/sd^3;	% original skewness
snrk = snr(sk, sk-s); 
sk = s*(1-p) + sk*p;

% Define the coefficients of the numerator (A*lam^3+B*lam^2+C*lam+D)

A=m(6)-3*sd*s*m(5)+3*sd^2*(s^2-1)*m(4)+sd^6*(2+3*s^2-s^4);
B=3*(m(5)-2*sd*s*m(4)+sd^5*s^3);
C=3*(m(4)-sd^4*(1+s^2));
D=s*sd^3;

a(7)=A^2;
a(6)=2*A*B;
a(5)=B^2+2*A*C;
a(4)=2*(A*D+B*C);
a(3)=C^2+2*B*D;
a(2)=2*C*D;
a(1)=D^2;

% Define the coefficients of the denominator (A2+B2*lam^2)

A2=sd^2;
B2=m(4)-(1+s^2)*sd^4;

b=zeros(1,7);
b(7)=B2^3;
b(5)=3*A2*B2^2;
b(3)=3*A2^2*B2;
b(1)=A2^3;


if 0, % test

	lam = -2:0.02:2;
	S = (A*lam.^3+B*lam.^2+C*lam+D)./...
		sqrt(b(7)*lam.^6 + b(5)*lam.^4 + b(3)*lam.^2 + b(1));
%	grd = ch.^2 - m(2) - sd * s * ch;
%	for lam = -1:0.01:1,
%		n = lam*100+101;
%		chp = ch + lam*grd;
%		S2(n) = mean2(chp.^3)/abs(mean2(chp.^2))^(1.5);
%	end
	lam = -2:0.02:2;
figure(1);plot(lam,S);grid;drawnow
%	snr(S2, S-S2)

end % test

% Now I compute its derivative with respect to lambda

d(8) = B*b(7);
d(7) = 2*C*b(7) - A*b(5);
d(6) = 3*D*b(7);
d(5) = C*b(5) - 2*A*b(3);
d(4) = 2*D*b(5) - B*b(3);
d(3) = -3*A*b(1);
d(2) = D*b(3) - 2*B*b(1);
d(1) = -C*b(1);

d = d(8:-1:1);
mMlambda = roots(d);

tg = imag(mMlambda)./real(mMlambda);
mMlambda = real(mMlambda(find(abs(tg)<1e-6)));
lNeg = mMlambda(find(mMlambda<0));
if length(lNeg)==0,
	lNeg = -1/eps;
end
lPos = mMlambda(find(mMlambda>=0));
if length(lPos)==0,
        lPos = 1/eps;
end
lmi = max(lNeg);
lma = min(lPos);

lam = [lmi lma];
mMnewSt = polyval([A B C D],lam)./(polyval(b(7:-1:1),lam)).^0.5;
skmin = min(mMnewSt);
skmax = max(mMnewSt);


% Given a desired skewness, solves for lambda

if sk<=skmin & Warn,
        lam = lmi;
        warning('Saturating (down) skewness!');
        skmin
elseif sk>=skmax & Warn,
        lam = lma;
        warning('Saturating (up) skewness!');
        skmax
else


% The equation is sum(c.*lam.^(0:6))=0

c=a-b*sk^2;

c=c(7:-1:1);

r=roots(c);

% Chose the real solution with minimum absolute value with the rigth sign
lam=-Inf;
co=0;
for n=1:6,
	tg = imag(r(n))/real(r(n));
	if (abs(tg)<1e-6)&(sign(real(r(n)))==sign(sk-s)),
		co=co+1;
		lam(co)=real(r(n));
	end
end
if min(abs(lam))==Inf & Warn,
	display('Warning: Skew adjustment skipped!');
	lam=0;
end

p=[A B C D];

if length(lam)>1,
	foo=sign(polyval(p,lam));
	if any(foo==0),
		lam = lam(find(foo==0));
	else
		lam = lam(find(foo==sign(sk)));		% rejects the symmetric solution
	end
	if length(lam)>0,
		lam=lam(find(abs(lam)==min(abs(lam))));	% the smallest that fix the skew
		lam=lam(1);
	else
		lam = 0;
	end
end
end % if else

% Modify the channel
chm=ch+lam*(ch.^2-sd^2-sd*s*ch);	% adjust the skewness
chm=chm*sqrt(m(2)/mean2(chm.^2));		% adjust the variance
chm=chm+me;				% adjust the mean
					% (These don't affect the skewness)
% Check the result
%mem=mean2(chm);
%sk2=mean2((chm-mem).^3)/mean2((chm-mem).^2).^(3/2);
%sk - sk2
%SNR=snr(sk,sk-sk2)




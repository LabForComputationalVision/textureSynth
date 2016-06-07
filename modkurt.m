function [chm, snrk] = modkurt(ch,k,p);

% Modify the kurtosis in one step, by moving in gradient direction until
% reaching the desired kurtosis value. 
% It does not affect the mean nor the variance, but it affects the skewness.
% This operation is not an orthogonal projection, but the projection angle is
% near pi/2 when k is close to the original kurtosis, which is a realistic assumption
% when doing iterative projections in a pyramid, for example (small corrections
% to the channels' statistics). 
%
% [chm, snrk] = modkurt(ch,k,p);
%	ch: channel
%	k: desired kurtosis (k=M4/M2^2)	
%       p [OPTIONAL]:   mixing proportion between k0 and k
%                       it imposes (1-p)*k0 + p*k,
%                       being k0 the current kurtosis.
%                       DEFAULT: p = 1;


% Javier Portilla,  Oct.12/97, NYU

Warn = 0;  % Set to 1 if you want to see warning messages
if ~exist('p'),
  p = 1;
end

me=mean2(ch);
ch=ch-me;

% Compute the moments

m=zeros(12,1);
for n=2:12,
	m(n)=mean2(ch.^n);
end

% The original kurtosis

k0=m(4)/m(2)^2;
snrk = snr(k, k-k0);
if snrk > 60,
	chm = ch+me;
	return
end
k = k0*(1-p) + k*p;

% Some auxiliar variables

a=m(4)/m(2);

% Coeficients of the numerator (A*lam^4+B*lam^3+C*lam^2+D*lam+E)

A=m(12)-4*a*m(10)-4*m(3)*m(9)+6*a^2*m(8)+12*a*m(3)*m(7)+6*m(3)^2*m(6)-...
	4*a^3*m(6)-12*a^2*m(3)*m(5)+a^4*m(4)-12*a*m(3)^2*m(4)+...
	4*a^3*m(3)^2+6*a^2*m(3)^2*m(2)-3*m(3)^4;
B=4*(m(10)-3*a*m(8)-3*m(3)*m(7)+3*a^2*m(6)+6*a*m(3)*m(5)+3*m(3)^2*m(4)-...
	a^3*m(4)-3*a^2*m(3)^2-3*m(4)*m(3)^2);
C=6*(m(8)-2*a*m(6)-2*m(3)*m(5)+a^2*m(4)+2*a*m(3)^2+m(3)^2*m(2));
D=4*(m(6)-a^2*m(2)-m(3)^2);
E=m(4);

% Define the coefficients of the denominator (F*lam^2+G)^2

F=D/4;
G=m(2);

% test
test = 0;

if test,

        grd = ch.^3 - a*ch - m(3);
        lam = -0.001:0.00001:0.001;
        k = (A*lam.^4+B*lam.^3+C*lam.^2+D*lam+E)./...
                (F*lam.^2 + G).^2;
        for lam = -0.001:0.00001:0.001,
                n = lam*100000+101;
                chp = ch + lam*grd;
                k2(n) = mean2(chp.^4)/mean2(chp.^2)^2;
                %k2(n) = mean2(chp.^4);
        end
        lam = -0.001:0.00001:0.001;
        snr(k2, k-k2)

end % test

% Now I compute its derivative with respect to lambda
% (only the roots of derivative = 0 )

d(1) = B*F;
d(2) = 2*C*F - 4*A*G;
d(3) = 4*F*D -3*B*G - D*F;
d(4) = 4*F*E - 2*C*G;
d(5) = -D*G;

mMlambda = roots(d);

tg = imag(mMlambda)./real(mMlambda);
mMlambda = mMlambda(find(abs(tg)<1e-6));
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
mMnewKt = polyval([A B C D E],lam)./(polyval([F 0 G],lam)).^2;
kmin = min(mMnewKt);
kmax = max(mMnewKt);

% Given a desired kurtosis, solves for lambda

if k<=kmin & Warn,
        lam = lmi;
        warning('Saturating (down) kurtosis!');
        kmin
elseif k>=kmax & Warn,
        lam = lma;
        warning('Saturating (up) kurtosis!');
        kmax
else

% Coeficients of the algebraic equation

c0 = E - k*G^2;
c1 = D;
c2 = C - 2*k*F*G;
c3 = B;
c4 = A - k*F^2;

% Solves the equation

r=roots([c4 c3 c2 c1 c0]);

% Chose the real solution with minimum absolute value with the rigth sign

tg = imag(r)./real(r);
%lambda = real(r(find(abs(tg)<1e-6)));
lambda = real(r(find(abs(tg)==0)));
if length(lambda)>0,
	lam = lambda(find(abs(lambda)==min(abs(lambda))));
	lam = lam(1);
else
	lam = 0;
end

end % if ... else


% Modify the channel

chm=ch+lam*(ch.^3-a*ch-m(3));	% adjust the kurtosis
chm=chm*sqrt(m(2)/mean2(chm.^2));	% adjust the variance
chm=chm+me;				% adjust the mean

% Check the result
%k2=mean2((chm-me).^4)/(mean2((chm-me).^2))^2;
%SNR=snr(k,k-k2)








 

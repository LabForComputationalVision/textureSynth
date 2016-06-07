% [newX, snr1, snr2, Mx, My] = adjustCorr2s(X, Cx, Y, Cxy, MODE, p) 
%
% Linearly adjust variables in X to have correlation Cx, and cross-correlation Cxy.
% Rows of X, Y, and newX are samples of (random) row-vectors, such that:
%   1:  newX = X * Mx + Y * My
%   2:  newX' * newX = Cx
%   3:  newX' * Y = Cxy
%
% MODE is optional:
%   0 => choose randomly from the space of linear solutions
%   1 => simplest soln
%   2 => minimize angle change
%   3 => Simple rotational (DEFAULT) 
%   4 => SVD minimal vector change soln
%
% p is optional:
%   Imposes an intermediate value of correlation between the current ones
%   Bx and Bxy and the specified Cx and Cxy:
%	Cx' = (1-p)*Bx + p*Cx;
%	Cxy' = (1-p)*Bxy + p*Cxy;
%   DEFAULT is p=1.


% EPS, 11/25/97

function [newX,snr1,snr2,Mx,My] = adjustCorr2s(X, Cx, Y, Cxy, mode, p)

Warn = 0; % Set to 1 if you want to display warning messages
if (exist('mode') ~= 1)
  mode = 3;
end
if (exist('p') ~= 1)
  p = 1;
end

Bx = innerProd(X) / size(X,1);
Bxy = (X' * Y) / size(X,1);
By = innerProd(Y) / size(X,1);
iBy = inv(By);

Current = Bx - (Bxy * iBy * Bxy');
Cx0 = Cx;
Cx = (1-p)*Bx + p*Cx;
Cxy0 = Cxy;
Cxy = (1-p)*Bxy + p*Cxy;
Desired = Cx - (Cxy * iBy * Cxy');

[E, D] = eig(Current);
D = diag(D);
if any(D < 0) & Warn
  ind = find(D<0);
  fprintf(1,'Warning: negative current eigenvalues: %d\n',D(ind)');
end
[junk,Ind] = sort(D);
D = diag(sqrt(D(Ind(size(Ind,1):-1:1))));
E = E(:,Ind(size(Ind,1):-1:1));

[Eo,Do] = eig(Desired);
Do = diag(Do);
if any(Do < 0) & Warn
  ind = find(Do<0);
  fprintf(1,'Warning: negative desired eigenvalues: %d\n',Do(ind)');
end
[junk,Ind] = sort(Do);
Do = diag(sqrt(Do(Ind(size(Ind,1):-1:1))));
Eo = Eo(:,Ind(size(Ind,1):-1:1));

if (mode == 0)
  Orth = orth(rand(size(D)));
elseif (mode == 1) % eye
  Orth = eye(size(D));
elseif (mode == 2) % simple
  A = [ eye(size(Cx)); -iBy*Bxy' ];
  Ao =  [ eye(size(Cx)); -iBy*Cxy' ];
  [U,S,V] = svd(E' * pinv(A) * Ao * Eo);
  Orth = U * V';
elseif (mode == 3)
  Orth = E' * Eo;
else     % SVD
  A = [ eye(size(Cx)); -iBy*Bxy' ];
  Ao =  [ eye(size(Cx)); -iBy*Cxy' ];
  [U,S,V] = svd(D * E' * pinv(A) * Ao * Eo * inv(Do));
  Orth = U * V';
end

Mx =  E * inv(D) * Orth * Do * Eo';
My =  iBy * (Cxy' - Bxy' * Mx);
newX = X * Mx + Y * My;

if Cx0~=Bx,
	snr1=10*log10(sum(sum(Cx0.^2))/sum(sum((Cx0-Bx).^2)));
else
	snr1 = Inf;
end
if Cxy0~=Bxy,
	snr2=10*log10(sum(sum(Cxy0.^2))/sum(sum((Cxy0-Bxy).^2)));
else
	snr2 = Inf;
end

% [newX, snr1, M] = adjustCorr1s(X, Cx, MODE, p) 
%
% Linearly adjust variables in X to have correlation Cx.
% Rows of X and newX are samples of a (random) row-vector, such that:
%    1:  newX = X * M    
%    2:  newX' * newX = Cx 
%
% MODE is optional:
%   0 => choose randomly from the space of linear solutions
%   1 => simplest soln
%   2 => minimize angle change (DEFAULT) 
%   3 => SVD minimal vector change soln
%
% p is optional:
%   Imposes an intermediate value of correlation between the current one
%   C and Cx:
%	Cx' = (1-p)*C + p*Cx;
%   DEFAULT is p=1.

%  EPS, 11/23/97.

function [newX, snr1, M] = adjustCorr1s(X,Co,mode,p)

if (exist('mode') ~= 1)
  mode = 2;
end

if (exist('p') ~= 1)
  p = 1;
end

C = innerProd(X) / size(X,1);
[E, D] = eig(C);
D = diag(D);
[junk,Ind] = sort(D);
D = diag(sqrt(D(Ind(size(Ind,1):-1:1))));
E = E(:,Ind(size(Ind,1):-1:1));

Co0 = Co;
Co = (1-p)*C + p*Co;

[Eo,Do] = eig(Co);
Do = diag(Do);
[junk,Ind] = sort(Do);
Do = diag(sqrt(Do(Ind(size(Ind,1):-1:1))));
Eo = Eo(:,Ind(size(Ind,1):-1:1));

if (mode == 0)
  Orth = orth(rand(size(C)));
elseif (mode == 1) % eye
  Orth = eye(size(C));
elseif (mode == 2) % simple
  Orth = E' * Eo;
else     % SVD
  [U,S,V] = svd(D * E' * Eo * inv(Do));
  Orth = U * V';
end

M =  E * inv(D) * Orth * Do * Eo';

newX = X * M;

snr1=10*log10(sum(sum(Co0.^2))/sum(sum((Co0-C).^2)));

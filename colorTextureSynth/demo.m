
im0 = imread('olives256.o.bmp');
im0 = im0(101+15:228+15,101-28:228-28,:);

Nsx = 128;  % Synthetic image dimensions
Nsy = 128;

Nsc = 3; % Number of pyramid scales
Nor = 4; % Number of orientations
Na = 7; % Number of spatial neighbors considered for spatial correlations
Niter = 25; % Number of iterations of the synthesis loop

[params] = textureColorAnalysis(im0, Nsc, Nor, Na);
tic; im = textureColorSynthesis(params, [Nsy Nsx], Niter); toc


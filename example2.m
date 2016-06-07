% Example 2:  Seamless blending of real and synthetic texture in an
% image, using Portilla-Simoncelli texture analysis/synthesis code,
% based on alternate projections onto statistical constraints in a
% complex overcomplete wavelet representation.
%
% See Readme.txt, and headers of textureAnalysis.m and
% textureSynthesis.m for more details.
%
% Javier Portilla (javier@decsai.ugr.es).  March, 2001

close all

Nsc = 4; % Number of scales
Nor = 4; % Number of orientations
Na = 5;  % Spatial neighborhood is Na x Na coefficients
			% It must be an odd number!

im0 = pgmRead('nuts.pgm');	% Warning: im0 is a double float matrix!

params = textureAnalysis(im0, Nsc, Nor, Na);

Niter = 25;	% Number of iterations of synthesis loop
Nsx = 192;	% Size of synthetic image is Nsy x Nsx
Nsy = 192;	% Warning: both dimensions must be multiple of 2^(Nsc+2)

% Use a mask and the original image to synthesize an image with the
% left side synthetic and the right side real data.
% The effective mask is M = (mask>0), its smoothness is for avoiding
% border effects.
ramp = meshgrid(1:Nsx/4,1:Nsy)*4/Nsy;
mask = [zeros(Nsy,Nsx/2) ramp ramp(:,Nsx/4:-1:1)];
mask =  1/2*(1-cos(pi*mask));

showIm(double(mask>0), 'auto', 'auto', 'Mask');
display('Press any key to continue');pause

imKeep = zeros(Nsx*Nsy,2);
imKeep(:,1) = reshape(mask, [Nsy*Nsx,1]);
imKeep(:,2) = reshape(im0(1:Nsy,1:Nsx), [Nsy*Nsx,1]);	% Original

res = textureSynthesis(params, [Nsy Nsx], Niter,[],imKeep);

close all
figure(1);showIm(double(mask>0), 'auto', 'auto', 'Mask');
figure(2);showIm(im0, 'auto', 'auto', 'Original Texture');
figure(3);showIm(res, 'auto', 'auto', 'Blended Original and Synthetic Texture');

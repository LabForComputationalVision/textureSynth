function [im,snrP,imS] = textureSynthesis(params, im0, Niter, cmask, imask)

% [res,snrP,imS] = textureSynthesis(params, initialIm, Niter, cmask, imask)
%
% Synthesize texture applying Portilla-Simoncelli model/algorithm.
%	
% params: structure containing texture parameters (as returned by textureAnalysis).
%
% im0: initial image, OR a vector (Ydim, Xdim, [SEED]) containing
% dimensions of desired image and an optional seed for the random
% number generator.  If dimensions are passed, initial image is
% Gaussian white noise.
%
% Niter (optional): Number of iterations.  Default = 50.
%
% cmask (optional): binary column vector (4x1) indicating which sets of
% constraints we want to apply in the synthesis. The four sets are:
%               1) Marginal statistics (mean, var, skew, kurt, range)
%               2) Correlation of subbands (space, orientation, scale)
%               3) Correlation of magnitude responses (sp, or, sc)
%               4) Relative local phase
%
% imask (optional): imsizex2 matrix.  First column is a mask, second
% 	column contains the image values to be imposed. If only one column is
%	provided, it assumes it corresponds to the image values, and it uses
%	a raised cosine square for the mask.
% snrP (optional):	Set of adjustment values (in dB) of the parameters.
% imS (optional):	Sequence of synthetic images, from niter = 1 to 2^n, being
% 			n = floor(log2(Niter)).

% Javier Portilla and Eero Simoncelli.
% Work described in:
%  "A Parametric Texture Model based on Joint Statistics of Complex Wavelet Coefficients".
%  J Portilla and E P Simoncelli. Int'l Journal of Computer Vision,
%  vol.40(1), pp. 49-71, Dec 2000.   
%
% Please refer to this publication if you use the program for research or
% for technical applications. Thank you.
%
% Copyright, Center for Neural Science, New York University, January 2001.
% All rights reserved.

Warn = 0;  % Set to 1 if you want to see warning messages

%% Check required args are passed:
if (nargin < 2)
  error('Function called with too few input arguments');
end

if ( ~exist('Niter') | isempty(Niter) )
  Niter = 50;
end

if (exist('cmask') & ~isempty(cmask) )
  cmask = (cmask > 0.5);  % indices of ones in mask
else
  cmask = ones(4,1);
end

%% Extract parameters  
statg0 = params.pixelStats;
mean0 = statg0(1); var0 =  statg0(2);
skew0 =  statg0(3); kurt0 =  statg0(4);
mn0 =  statg0(5);  mx0 = statg0(6);
statsLPim = params.pixelLPStats;
skew0p = statsLPim(:,1);
kurt0p = statsLPim(:,2);
vHPR0 = params.varianceHPR;
acr0 = params.autoCorrReal;
ace0 = params.autoCorrMag;
magMeans0 = params.magMeans;
C0 = params.cousinMagCorr;
Cx0 = params.parentMagCorr;
Crx0 = params.parentRealCorr;

%% Extract {Nsc, Nor, Na} from params
tmp = size(params.autoCorrMag);
Na = tmp(1); Nsc = tmp(3);
Nor = tmp(length(tmp))*(length(tmp)==4) + (length(tmp)<4);
la = (Na-1)/2;

%% If im0 is a vector of length 2, create Gaussian white noise image of this
%% size, with desired pixel mean and variance.  If vector length is
%% 3,  use the 3rd element to seed the random number generator.
if ( length(im0) <= 3 )
  if ( length(im0) == 3)
    randn('state', im0(3)); % Reset Seed
    im0 = im0(1:2);
  end
  im = mean0 + sqrt(var0)*randn(im0);
else
  im = im0;
end

%% If the spatial neighborhood Na is too big for the lower scales,
%% "modacor22.m" will make it as big as the spatial support at
%% each scale:
[Ny,Nx] = size(im);
nth = log2(min(Ny,Nx)/Na);
if nth<Nsc+1 & Warn,
  fprintf(1,'Warning: Na will be cut off for levels above #%d !\n',floor(nth));
end

if  exist('imask') & ~isempty(imask),
	if  size(imask,1) ~= prod(size(im))
  		error(sprintf('imask size %d does not match image dimensions [%d,%d]',...
		size(imask,1), size(im,1), size(im,2)));
	end
	if size(imask,2) == 1,
        	mask = (cos(-pi/2:2*pi/Ny:pi*(1-2/Ny)/2)).'*cos(-pi/2:2*pi/Nx:pi*(1-2/Nx)/2);
        	mask = mask.^2;
        	aux = zeros(size(im));
        	aux(Ny/4+1:Ny/4+Ny/2,Nx/4+1:Nx/4+Nx/2)=mask;
        	mask = aux;
	else
        	mask = reshape(imask(:,1),size(im));
	end
end

current_gcf = gcf;

imf = max(1,current_gcf.Number - 1); snrf = imf+1;
figure(imf);  clf
subplot(1,2,1); grayRange = showIm(im,'auto',1); title('Starting image');
drawnow

prev_im=im;

current_gcf = gcf;

imf = max(1, current_gcf.Number - 1);
figure(imf);   
clf;showIm(im,'auto',1); title(sprintf('iteration 0'));

nq = 0;
Nq = floor(log2(Niter));
imS = zeros(Ny,Nx,Nq);

%% MAIN LOOP
for niter = 1:Niter

%p = niter/Niter; 
p = 1;

  %% Build the steerable pyramid
  [pyr,pind] = buildSCFpyr(im,Nsc,Nor-1);

  if ( any(vectify(mod(pind,4))) )
    error('Algorithm will fail: band dimensions are not all multiples of 4!');
  end

  %% Subtract mean of lowBand:
  nband = size(pind,1);
  pyr(pyrBandIndices(pind,nband)) = ...
      pyrBand(pyr,pind,nband) - mean2(pyrBand(pyr,pind,nband));

  apyr = abs(pyr);

  %% Adjust autoCorr of lowBand
  nband = size(pind,1); 
  ch = pyrBand(pyr,pind,nband);
  Sch = min(size(ch)/2);
  nz = sum(sum(~isnan(acr0(:,:,Nsc+1))));
  lz = (sqrt(nz)-1)/2;
  le = min(Sch/2-1,lz);
  im = real(ch);  %Reconstructed image: initialize to lowband
  [mpyr,mpind] = buildSFpyr(im,0,0);
  im = pyrBand(mpyr,mpind,2);
  vari =  acr0(la+1:la+1,la+1:la+1,Nsc+1);
if cmask(2),
  if vari/var0 > 1e-4,
  	[im, snr2(niter,Nsc+1)] = ...
      	modacor22(im, acr0(la-le+1:la+le+1,la-le+1:la+le+1,Nsc+1),p);
  else
	im = im*sqrt(vari/var2(im));
  end
  if (var2(imag(ch))/var2(real(ch)) > 1e-6)
    fprintf(1,'Discarding non-trivial imaginary part, lowPass autoCorr!');
  end
  im = real(im);
end % cmask(2)
if cmask(1),
  if vari/var0 > 1e-4,
  	[im,snr7(niter,2*(Nsc+1)-1)] = modskew(im,skew0p(Nsc+1),p);	% Adjusts skewness 
  	[im,snr7(niter,2*(Nsc+1))] = modkurt(im,kurt0p(Nsc+1),p);	% Adjusts kurtosis 
  end
end	% cmask(2)
 
  %% Subtract mean of magnitude
if cmask(3),
  magMeans = zeros(size(pind,1), 1);
  for nband = 1:size(pind,1)
    indices = pyrBandIndices(pind,nband);
    magMeans(nband) = mean2(apyr(indices));
    apyr(indices) = apyr(indices) - magMeans(nband);
  end
end	% cmask(3)

  %% Coarse-to-fine loop:
  for nsc = Nsc:-1:1
    
    firstBnum = (nsc-1)*Nor+2;
    cousinSz = prod(pind(firstBnum,:));
    ind = pyrBandIndices(pind,firstBnum);
    cousinInd = ind(1) + [0:Nor*cousinSz-1];

    %% Interpolate parents
if (cmask(3) | cmask(4)),
    if (nsc<Nsc)
      parents = zeros(cousinSz,Nor);
      rparents = zeros(cousinSz,Nor*2);
      for nor = 1:Nor
	nband = (nsc+1-1)*Nor+nor+1; 

        tmp = expand(pyrBand(pyr, pind, nband),2)/4;
	rtmp = real(tmp); itmp = imag(tmp);
        tmp = sqrt(rtmp.^2 + itmp.^2) .* exp(2 * sqrt(-1) * atan2(rtmp,itmp));
        rparents(:,nor) = vectify(real(tmp));
        rparents(:,Nor+nor) = vectify(imag(tmp));

        tmp = abs(tmp);
        parents(:,nor) = vectify(tmp - mean2(tmp));
      end
    else
      rparents = [];
      parents = [];
    end
end % if (cmask(3) | cmask(4))

if cmask(3),
    %% Adjust cross-correlation with MAGNITUDES at other orientations/scales:
    cousins = reshape(apyr(cousinInd), [cousinSz Nor]);
    nc = size(cousins,2);   np = size(parents,2);
    if (np == 0)    
      [cousins, snr3(niter,nsc)] = adjustCorr1s(cousins, C0(1:nc,1:nc,nsc), 2, p);
    else
      [cousins, snr3(niter,nsc), snr4(niter,nsc)] = ...
	  adjustCorr2s(cousins, C0(1:nc,1:nc,nsc), parents, Cx0(1:nc,1:np,nsc), 3, p);
    end
    if (var2(imag(cousins))/var2(real(cousins)) > 1e-6)
      fprintf(1,'Non-trivial imaginary part, mag crossCorr, lev=%d!\n',nsc);
    else
      cousins = real(cousins);
      ind = cousinInd;
      apyr(ind) = vectify(cousins);
    end

    %% Adjust autoCorr of mag responses 
    nband = (nsc-1)*Nor+2; 
    Sch = min(pind(nband,:)/2);
    nz = sum(sum(~isnan(ace0(:,:,nsc,1))));
    lz = (sqrt(nz)-1)/2;
    le = min(Sch/2-1,lz);
    for nor = 1:Nor,
      nband = (nsc-1)*Nor+nor+1; 
      ch = pyrBand(apyr,pind,nband);
      [ch, snr1(niter,nband-1)] = modacor22(ch,...
	ace0(la-le+1:la+le+1,la-le+1:la+le+1,nsc,nor), p);
      ch = real(ch);
      ind = pyrBandIndices(pind,nband);
      apyr(ind) = ch;
      %% Impose magnitude:
      mag = apyr(ind) + magMeans0(nband);
      mag = mag .* (mag>0);
      pyr(ind) = pyr(ind) .* (mag./(abs(pyr(ind))+(abs(pyr(ind))<eps)));
    end
end   % cmask(3)

    %% Adjust cross-correlation of REAL PARTS at other orientations/scales:
    cousins = reshape(real(pyr(cousinInd)), [cousinSz Nor]);
    Nrc = size(cousins,2);  Nrp = size(rparents,2);
    if  cmask(4) & (Nrp ~= 0)
	a3 = 0; a4 = 0;
        for nrc = 1:Nrc,
          cou = cousins(:,nrc);
          [cou, s3, s4] = ...
            adjustCorr2s(cou,mean(cou.^2),rparents,Crx0(nrc,1:Nrp,nsc), 3);
	  a3 = s3 + a3;
	  a4 = s4 + a4;
          cousins(:,nrc) = cou;
        end
	snr4r(niter,nsc) = a4/Nrc;
    end
    if (var2(imag(cousins))/var2(real(cousins)) > 1e-6)
      fprintf(1,'Non-trivial imaginary part, real crossCorr, lev=%d!\n',nsc);
    else
      %%% NOTE: THIS SETS REAL PART ONLY - signal is now NONANALYTIC!
      pyr(cousinInd) = vectify(cousins(1:Nor*cousinSz));
    end

    %% Re-create analytic subbands
    dims = pind(firstBnum,:);
    ctr = ceil((dims+0.5)/2);
    ang = mkAngle(dims, 0, ctr);
    ang(ctr(1),ctr(2)) = -pi/2;
    for nor = 1:Nor,
      nband = (nsc-1)*Nor+nor+1; 
      ind = pyrBandIndices(pind,nband); 
      ch = pyrBand(pyr, pind, nband);
      ang0 = pi*(nor-1)/Nor;
      xang = mod(ang-ang0+pi, 2*pi) - pi;
      amask = 2*(abs(xang) < pi/2) + (abs(xang) == pi/2);
      amask(ctr(1),ctr(2)) = 1;
      amask(:,1) = 1;
      amask(1,:) = 1; 
      amask = fftshift(amask);
      ch = ifft2(amask.*fft2(ch));	% "Analytic" version
      pyr(ind) = ch;
    end

    %% Combine ori bands
    bandNums = [1:Nor] + (nsc-1)*Nor+1;  %ori bands only
    ind1 = pyrBandIndices(pind, bandNums(1));
    indN = pyrBandIndices(pind, bandNums(Nor));
    bandInds = [ind1(1):indN(length(indN))];
    %% Make fake pyramid, containing dummy hi, ori, lo
    fakePind = pind([bandNums(1), bandNums, bandNums(Nor)+1],:);
    fakePyr = [zeros(prod(fakePind(1,:)),1);...
	 real(pyr(bandInds)); zeros(prod(fakePind(size(fakePind,1),:)),1)]; 
    ch = reconSFpyr(fakePyr, fakePind, [1]);     % recon ori bands only
    im = real(expand(im,2))/4;
    im = im + ch;
    vari =  acr0(la+1:la+1,la+1:la+1,nsc);
if cmask(2),
    if vari/var0 > 1e-4,
  	[im, snr2(niter,nsc)] = ...
      	modacor22(im, acr0(la-le+1:la+le+1,la-le+1:la+le+1,nsc), p);
    else
	im = im*sqrt(vari/var2(im));
    end
end	% cmask(2)
    im = real(im);  

if cmask(1),
  %% Fix marginal stats 
  if vari/var0 > 1e-4,
        [im,snr7(niter,2*nsc-1)] = modskew(im,skew0p(nsc),p);       % Adjusts skewness 
        [im,snr7(niter,2*nsc)] = modkurt(im,kurt0p(nsc),p);         % Adjusts kurtosis
  end
end	% cmask(1)

  end  %END Coarse-to-fine loop

  %% Adjust variance in HP, if higher than desired
if (cmask(2)|cmask(3)|cmask(4)),
  ind = pyrBandIndices(pind,1);
  ch = pyr(ind);
  vHPR = mean2(ch.^2);
  if vHPR > vHPR0,
	ch = ch * sqrt(vHPR0/vHPR);
	pyr(ind) = ch;
  end
end % cmask
  im = im + reconSFpyr(real(pyr), pind, [0]);  %recon hi only

  %% Pixel statistics
  means = mean2(im);
  vars = var2(im, means);
  snr7(niter,2*(Nsc+1)+1) = snr(var0,var0-vars);
  im = im-means;  			% Adjusts mean and variance
  [mns mxs] = range2(im + mean0);
  snr7(niter,2*(Nsc+1)+2) = snr(mx0-mn0,sqrt((mx0-mxs)^2+(mn0-mns)^2));
if cmask(1),
  im = im*sqrt(((1-p)*vars + p*var0)/vars);	
end	% cmaks(1)
  im = im+mean0;
if cmask(1),
  [im, snr7(niter,2*(Nsc+1)+3)] = modskew(im,skew0,p);	% Adjusts skewness (keep mean and variance)
  [im, snr7(niter,2*(Nsc+1)+4)] = modkurt(im,kurt0,p);	% Adjusts kurtosis (keep mean and variance,
					% but not skewness)
  im = max(min(im,(1-p)*max(max(im))+p*mx0),...
		  (1-p)*min(min(im))+p*mn0);		% Adjusts range (affects everything)	
else 
  snr7(niter,2*(Nsc+1)+3) = snr(skew0,skew0-skew2(im));
  snr7(niter,2*(Nsc+1)+4) = snr(kurt0,kurt0-kurt2(im));
end	% cmask(1)

  %% Force pixels specified by image mask
  if (exist('imask') & ~isempty(imask) )
    	im = mask.*reshape(imask(:,2 - (size(imask,2)==1)),size(im)) + ...
	 	(1-mask).*im;
  end

  snr6(niter,1) = snr(im-mean0,im-prev_im);

  if floor(log2(niter))==log2(niter),
	nq = nq + 1;
  	imS(:,:,nq) = im;
  end

  tmp = prev_im;
  prev_im=im;	

  figure(imf);
  subplot(1,2,1);
  showIm(im-tmp,'auto',1); title('Change');
  subplot(1,2,2);
  showIm(im,'auto',1); title(sprintf('iteration %d/%d',niter,Niter));
  drawnow
  
  % accelerator
  alpha = 0.8;
  im = im + alpha*(im - tmp);
  
commented = 1;  % set it to 0 for displaying convergence of parameters in SNR (dB)
if ~commented,    
   
% The graphs that appear reflect
% the relative distance of each parameter or group
% of parametersi, to the original's, in decibels.
% Note, however, that when the original parameters
% are close to zero, this measurement is meaningless.
% This is why in some cases it seems that some of
% the parameters do not converge at all.

figure(snrf);
if cmask(1)
  subplot(171); plot(snr7); title('Mrgl stats');
end
if cmask(2),
  subplot(172); plot(snr2); title('Raw auto');
end
if cmask(3),
  subplot(173); plot(snr1); title('Mag auto'); 
  subplot(174); plot(snr3); title('Mag ori');
  subplot(175); plot(snr4); title('Mag scale');
end
if (Nrp > 0) & cmask(4),
  subplot(176); plot(snr4r); title('Phs scale');
end
  subplot(177); plot(snr6); title('Im change');
  drawnow
  
end  % if ~commented

end %END  MAIN LOOP

im = prev_im;

snrP = [snr7 snr2 snr1 snr3 snr4 snr4r snr6];

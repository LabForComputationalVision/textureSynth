function [im,snrP,imS] = textureColorSynthesis(params, im0, Niter, cmask, imask, pmask)

% res = textureColorSynthesis(params, initialIm, Niter, cmask, imask, pmask)
%
% Synthesize color texture with Portilla-Simoncelli model/algorithm.   
%
% params: structure containing texture parameters (as returned by textureColorAnalysis).
%
% im0: initial image (true color, uint8), OR a vector (Ydim, Xdim, [SEED])
%       containing dimensions of desired image and an optional seed for the random
%       number generator.  If dimensions are passed, initial image is
%       Gaussian white noise.
%
% Niter (optional): Number of iterations.  Default = 50.
%
% cmask (optional): binary column vector (5x1) indicating which sets of
%       constraints we want to apply in the synthesis and the display option.
%       The four constraints are:
%               1) Marginal statistics (moments at different scales) 
%               2) Auto-Correlation of image (at different scales) 
%               3) Correlation of magnitude responses (space, orientation, scale)
%               4) Relative local phase
%       The fifth element is:
%               5) Display convergence plots while processing (1 = on/ 0 = off, default)
%
% imask (optional): imsize(:) x 4 matrix.  First column is a binary mask, second,
%       third and fourth columns contain the image values to be imposed (R,G,B).
%
% pmask (optional): pyrsize x 2 matrix.  First column is a binary mask, second
%       column contains the (real) pyramid coefficients to be imposed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See readme.txt for further information about installing/using this code.
% See copyright.txt for restrictions on usage.
%
% J. Portilla and E. P. Simoncelli
% portilla@io.cfmac.csic.es, eero.simoncelli@nyu.edu
%
% October 2000, New York University, New York
% Released Version 1.0: January 2009, CSIC, Madrid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nclr = 3;   % Number of colors

%% Check required args are passed:
if (nargin < 2)
  error('Function called with too few input arguments');
end

if ~exist('Niter') | isempty(Niter) 
  Niter = 50;
end

if exist('cmask') & ~isempty(cmask) 
  cmask = (cmask > 0.5);  % indices of ones in mask
else 
  cmask = ones(5,1);
  cmask(5) = 1;
end

if exist('imask') & ~isempty(imask) 
  imaskInd = find( imask(:,1) > 0.5 );  % indices of ones in mask
else
  imaskInd = [];
end

if exist('pmask')  & ~isempty(pmask) 
  pmaskInd = find( pmask(:,1) > 0.5 );  % indices of ones in mask
else
  pmaskInd = [];
end

%% Extract parameters:
statg0 = params.pixelStats.';
mean0 = statg0(1,:); var0 =  statg0(2,:);
skew0 =  statg0(3,:); kurt0 =  statg0(4,:);
mn0 =  statg0(5,:);  mx0 = statg0(6,:);
statgP = params.pixelStatsPCA.';
skewP =  statgP(1,:); kurtP =  statgP(2,:);
mnP =  statgP(3,:);  mxP = statgP(4,:);
statsLPim = params.pixelLPStats;
skew0p = statsLPim(1:size(statsLPim,1)/2,:);
kurt0p = statsLPim(size(statsLPim,1)/2+1:size(statsLPim,1),:);
acr0 = params.autoCorrReal;
ace0 = params.autoCorrMag;
magMeans0 = params.magMeans;
C0 = params.cousinMagCorr;
Cx0 = params.parentMagCorr;
Cr0 = params.cousinRealCorr;
Crx0 = params.parentRealCorr;
vHPR0 = params.varianceHPR;
Cclr0 = params.colorCorr;

Nclr = 3;

%% Extract {Nsc, Nor, Na} from params
tmp = size(params.autoCorrMag);
Na = tmp(1); Nsc = tmp(3); Nor = tmp(4);
la = (Na-1)/2;

%% If im0 is a vector of length 2, create Gaussian white noise image of this
%% size, with desired pixel mean and variance.  If vector length is
%% 3,  use the 3rd element to seed the random number generator.
if ( length(im0) <= 3 )
  if ( length(im0) == 3)
    randn('state', im0(3)); % Reset Seed
    im0 = im0(1:2);
  end
  Ny = im0(1);
  Nx = im0(2);
  im = zeros(Ny, Nx, Nclr);
  for clr = 1:Nclr,
        im(:,:,clr) = randn(Ny,Nx);
  end
  im = reshape(im, Ny*Nx, Nclr);
  im = adjustCorr1s(im, Cclr0);
  im = real(im);
  im = im + ones(Ny*Nx,Nclr)*diag(mean0);
  im = reshape(im, Ny, Nx, Nclr);
else
  im = double(im0);
end

if (( exist('imask') == 1 ) & ~isempty(imask) )
  imMask= double(reshape(imask(:,1),Ny,Nx));
  imOrig = double(reshape(imask(:,2:4),Ny,Nx,Nclr));
  for clr =1 : Nclr,
        im(:,:,clr) = imOrig(:,:,clr).*imMask + im(:,:,clr).*(1-imMask);
  end
end


%% If the spatial neighborhood Na is too big for the lower scales,
%% "modacor22.m" will make it as big as the spatial support at
%% each scale:
[Ny,Nx,Nclr] = size(im);
nth = log2(min(Ny,Nx)/Na);
if nth<Nsc+1,
  fprintf(1,'Warning: Na will be cut off for levels above #%d !\n',floor(nth));
end

if ( ~isempty(imaskInd) &  (size(imask,1) ~= prod(size(im))/3) )
  error(sprintf('imask size %d does not match image dimensions [%d,%d]',...
                size(imask,1), size(im,1), size(im,2)));
end

if cmask(5),  % display numerical fit to constraints
imf = max(1,gcf-1); snrf = imf+1;
figure(imf);  clf
subplot(1,2,1); showIm(im(:,:,1),[0 255],256/max(Ny,Nx));
image(uint8(max(min(im,255),0))); axis('image');axis('off');
title('Starting image');
drawnow
end % cmask(5)

[V,D] = eig(Cclr0);

prev_im = im;

nq = 0;
Nq = floor(log2(Niter));
imS = zeros(Ny,Nx,Nclr,Nq);

%% MAIN LOOP
for niter = 1:Niter

if ~cmask(5),   % Inform of number of iterations, when display option is off
	niter
end

  % Work with the PCA color components
  im = reshape(im,Ny*Nx,Nclr);
  im = im - ones(Ny*Nx,Nclr)*diag(mean(im));
  im = im*V*pinv(sqrt(D));
  im = reshape(im,Ny,Nx,Nclr);

  p = 1;

  %% Build the steerable pyramid
  for clr = 1:Nclr,
        [pyr(:,clr),pind] = buildSCFpyr(im(:,:,clr),Nsc,Nor-1);
  end

  if ( any(vectify(mod(pind,4))) )
    error('Algorithm will fail: band dimensions are not all multiples of 4!');
  end
  if ( ~isempty(pmaskInd) &  (size(pmask,1) ~= size(pyr,1)) )
    error(sprintf('pmask size %d does not match pyramid size %d',...
                  size(pmask,1), size(pyr,1)));
  end

  %% Subtract mean of lowBand:
  nband = size(pind,1);
  for clr = 1:Nclr,
    pyr(pyrBandIndices(pind,nband),clr) = ...
        real(pyr(pyrBandIndices(pind,nband),clr)) - ...
        mean(real(pyr(pyrBandIndices(pind,nband),clr)));
  end
  apyr = abs(pyr);

  %% Subtract mean of magnitude
if cmask(3),
  magMeans = zeros(size(pind,1), Nclr);
  for clr = 1:Nclr,
    for nband = 1:size(pind,1)
      indices = pyrBandIndices(pind,nband);
      magMeans(nband,clr) = mean(apyr(indices,clr));
      apyr(indices,clr) = apyr(indices,clr) - magMeans(nband,clr);
    end
  end
end % cmask(3)

  %% Adjust central autoCorr of lowBand and its skew and kurt
  %% Adjust (try) spatial-color correlation in the LPR (low pass residual)
  [Nly, Nlx] = size(pyrLow(pyr(:,1),pind));
  Nly = 2*Nly; Nlx = 2*Nlx;
  LPRcross = zeros(Nly*Nlx,5,Nclr);
  for clr = 1:Nclr,
        tmp = real(expand(pyrLow(pyr(:,clr),pind),2))/4;
        LPRcross(:,:,clr) = [vectify(tmp),...
              vectify(shift(tmp,[0 2])), vectify(shift(tmp,[0 -2])), ...
              vectify(shift(tmp,[2 0])), vectify(shift(tmp,[-2 0]))];
  end  % clr
  Nrp = size(LPRcross,2)*Nclr;
  LPRcross = reshape(LPRcross,[Nly*Nlx Nrp]);
  [LPRcross, snr3r(niter,Nsc+1)]=...
        adjustCorr1s(LPRcross,Cr0(1:Nrp,1:Nrp,Nsc+1), 2);
  LPRcross = reshape(LPRcross,[Nly*Nlx 5 Nclr]);
  indices = pyrBandIndices(pind,nband);
  for clr = 1:Nclr,
        LPRcross(:,:,clr) = [LPRcross(:,1,clr),...
                vectify(shift(reshape(LPRcross(:,2,clr),Nly,Nlx),[0 -2])),...
                vectify(shift(reshape(LPRcross(:,3,clr),Nly,Nlx),[0 2])), ...
                vectify(shift(reshape(LPRcross(:,4,clr),Nly,Nlx),[-2 0])),...
                vectify(shift(reshape(LPRcross(:,5,clr),Nly,Nlx),[2 0]))];
        aux = mean(LPRcross(:,:,clr),2);
        aux = reshape(aux,Nly,Nlx);
        pyr(indices,clr) = 4*vectify(real(shrink(aux,2)));
  end  % clr


  for clr = 1:Nclr,
    ch = pyrBand(pyr(:,clr),pind,nband);
    [Nly Nlx] = size(ch);
    Sch = min(size(ch)/2);
    nz = sum(sum(~isnan(acr0(:,:,Nsc+1,clr))));
    lz = (sqrt(nz)-1)/2;
    le = min(Sch/2-1,lz);
    [mpyr,mpind] = buildSFpyr(real(ch),0,0);
    im(1:Nly,1:Nlx,clr) = pyrBand(mpyr,mpind,2);
    vari =  acr0(la+1:la+1,la+1:la+1,Nsc+1,clr);
if cmask(2),
    if vari/D(clr,clr) > 1e-3,
        [im(1:Nly,1:Nlx,clr), snr2(niter,Nsc+1,clr)] = ...
          	modacor22(im(1:Nly,1:Nlx,clr),...
	   	acr0(la-le+1:la+le+1,la-le+1:la+le+1,Nsc+1,clr),p);
    else
        im(1:Nly,1:Nlx,clr) = im(1:Nly,1:Nlx,clr)*sqrt(vari/var2(im(1:Nly,1:Nlx,clr)));
    end
    if (var2(imag(ch))/var2(real(ch)) > 1e-3)
      fprintf(1,'Discarding non-trivial imaginary part, lowPass autoCorr!');
    end
    im(1:Nly,1:Nlx,clr) = real(im(1:Nly,1:Nlx,clr));
end % cmask(2)
if cmask(1),
    if vari/D(clr,clr) > 1e-3,
        [im(1:Nly,1:Nlx,clr),snr7(niter,2*(Nsc+1)-1,clr)] = ...
	  modskew(im(1:Nly,1:Nlx,clr),skew0p(Nsc+1,clr),p);     % Adjusts skewness
        [im(1:Nly,1:Nlx,clr),snr7(niter,2*(Nsc+1),clr)] = ...
	  modkurt(im(1:Nly,1:Nlx,clr),kurt0p(Nsc+1,clr),p);     % Adjusts kurtosis
    end
end     % cmask(1)
  end  % clr	
 
  %% Coarse-to-fine loop:
  for nsc = Nsc:-1:1

    firstBnum = (nsc-1)*Nor+2;
    cousinSz = prod(pind(firstBnum,:));
    ind = pyrBandIndices(pind,firstBnum);
    cousinInd = ind(1) + [0:Nor*cousinSz-1];

    cousins = zeros(cousinSz,Nor,Nclr);
    rcousins = zeros(cousinSz,Nor,Nclr);
    parents = zeros(cousinSz,Nor,Nclr);
    rparents = zeros(cousinSz,2*Nor,Nclr);

if (cmask(3) | cmask(4)),
    for clr = 1:Nclr,

    %% Interpolate parents
      if (nsc<Nsc)
        for nor = 1:Nor
          nband = (nsc+1-1)*Nor+nor+1;

          tmp = expand(pyrBand(pyr(:,clr), pind, nband),2)/4;
          rtmp = real(tmp); itmp = imag(tmp);
          tmp = sqrt(rtmp.^2 + itmp.^2) .* exp(2 * sqrt(-1) * atan2(rtmp,itmp));
          rparents(:,nor,clr) = vectify(real(tmp));
          rparents(:,Nor+nor,clr) = vectify(imag(tmp));

          tmp = abs(tmp);
          parents(:,nor,clr) = vectify(tmp - mean2(tmp));
        end
      else
        tmp = real(expand(pyrLow(pyr(:,clr),pind),2))/4;
        rparents(:,1:5,clr) = [vectify(tmp),...
                  vectify(shift(tmp,[0 2])), vectify(shift(tmp,[0 -2])), ...
                  vectify(shift(tmp,[2 0])), vectify(shift(tmp,[-2 0]))];
        parents = [];
      end
      cousins(:,:,clr) = reshape(apyr(cousinInd,clr), [cousinSz Nor]);
      rcousins(:,:,clr) = reshape(real(pyr(cousinInd,clr)), [cousinSz Nor]);


    end  % clr
end % if (cmask(3) | cmask(4))

if cmask(3),
    %% Adjust cross-correlation with MAGNITUDES at other orientations/scales:
    nc = prod(size(cousins))/cousinSz;
    np = prod(size(parents))/cousinSz;
    cousins = reshape(cousins,[cousinSz nc]);
    parents = reshape(parents,[cousinSz np]);
    if (np == 0)
      [cousins, snr3(niter,nsc)] = adjustCorr1s(cousins, C0(:,:,nsc),2);
    else
      [cousins, snr3(niter,nsc), snr4(niter,nsc)] = ...
         adjustCorr2s(cousins, C0(:,:,nsc),parents,Cx0(:,:,nsc),3);
    end
    cousins = reshape(cousins, cousinSz, Nor, Nclr);

    for clr = 1:Nclr,

      cou = cousins(:,:,clr);
      if (var2(imag(cou))/var2(real(cou)) > 1e-3)
        fprintf(1,'Non-trivial imaginary part, mag crossCorr, lev=%d!\n',nsc);
      else
        cou = real(cou);
        apyr(cousinInd,clr) = vectify(cou);
      end

      %% Adjust autoCorr of magnitude responses
      nband = (nsc-1)*Nor+2;
      Sch = min(pind(nband,:)/2);
      nz = sum(sum(~isnan(ace0(:,:,nsc,1,clr))));
      lz = (sqrt(nz)-1)/2;
      le = min(Sch/2-1,lz);
      for nor = 1:Nor,
        nband = (nsc-1)*Nor+nor+1;
        ch = pyrBand(apyr(:,clr),pind,nband);
        [ch, snr1(niter,nband-1,clr)] = modacor22(ch,...
          ace0(la-le+1:la+le+1,la-le+1:la+le+1,nsc,nor,clr));
        ch = real(ch);
        apyr(pyrBandIndices(pind,nband),clr) = vectify(ch);
        %% Impose magnitude:
        ind = pyrBandIndices(pind,nband);
        mag = apyr(ind,clr) + magMeans0(nband,clr);
        mag = mag .* (mag>0);
        pyr(ind,clr) = pyr(ind,clr)...
	 	 .* (mag./(abs(pyr(ind,clr))+(abs(pyr(ind,clr))<eps)));
      end
      rcousins(:,:,clr) = reshape(real(pyr(cousinInd,clr)), [cousinSz Nor]);

    end         % clr

end % if cmask(3)

    %% Adjust cross-correlation of REAL PARTS at other orientations/scales:
    if (nsc == Nsc),
        rparents = rparents(:,1:5,:);
    end
    Nrp = prod(size(rparents))/cousinSz;
    Nrc = prod(size(rcousins))/cousinSz;
    rparents = reshape(rparents,[cousinSz Nrp]);
    rcousins = reshape(rcousins,[cousinSz Nrc]);
 
    if ((Nrp == 0) & cmask(2))|(~cmask(4) & cmask(2)),
      [rcousins, snr3r(niter,nsc)]=...
        adjustCorr1s(rcousins,Cr0(1:Nrc,1:Nrc,nsc), 2);
    elseif (cmask(4) & cmask(2)),
        [rcousins, snr3r(niter,nsc),snr4r(niter,nsc)] = ...
                adjustCorr2s(rcousins,Cr0(1:Nrc,1:Nrc,nsc),...
                rparents,Crx0(1:Nrc,1:Nrp,nsc), 3);
    elseif (cmask(4) & ~cmask(2) & (nsc~=Nsc))
        for nrc=1:Nrc,
            cou = rcousins(:,nrc);
            [cou, snr3r(niter,nsc),snr4r(niter,nsc)] = ...
                adjustCorr2s(cou,Cr0(nrc,nrc,nsc),...
                rparents,Crx0(nrc,1:Nrp,nsc), 3);
            rcousins(:,nrc) = cou;
        end
    end
    rcousins = reshape(rcousins, cousinSz, Nor, Nclr);

    for clr = 1:Nclr,

if  cmask(4),
      if (var2(imag(rcousins(:,:,clr)))/var2(real(rcousins(:,:,clr))) > 1e-3)
        fprintf(1,'Non-trivial imaginary part, real crossCorr, lev=%d!\n',nsc);
       end
        %%% NOTE: THIS SETS REAL PART ONLY!
        pyr(cousinInd,clr) = vectify(real(rcousins(:,:,clr)));
end % cmask(4)

      %% Re-create analytic subbands
      dims = pind(firstBnum,:);
      ctr = ceil((dims+0.5)/2);
      ang = mkAngle(dims, 0, ctr);
      ang(ctr(1),ctr(2)) = -pi/2;
      for nor = 1:Nor,
        nband = (nsc-1)*Nor+nor+1;
        ind = pyrBandIndices(pind,nband);
        ch = pyrBand(pyr(:,clr), pind, nband);
        ang0 = pi*(nor-1)/Nor;
        xang = mod(ang-ang0+pi, 2*pi) - pi;
        amask = 2*(abs(xang) < pi/2) + (abs(xang) == pi/2);
        amask(ctr(1),ctr(2)) = 1;
        amask(:,1) = 1;
        amask(1,:) = 1;
        amask = fftshift(amask);
        ch = ifft2(amask.*fft2(ch));      % "Analytic" version
        pyr(ind,clr) = vectify(ch);
      end

      %% Combine ori bands
      bandNums = [1:Nor] + (nsc-1)*Nor+1;  %ori bands only
      ind1 = pyrBandIndices(pind, bandNums(1));
      indN = pyrBandIndices(pind, bandNums(Nor));
      bandInds = [ind1(1):indN(length(indN))];
      %% Make fake pyramid, containing dummy hi, ori, lo
      fakePind = [pind(bandNums(1),:);pind(bandNums(1):bandNums(Nor)+1,:)];
      fakePyr = [zeros(prod(fakePind(1,:)),1);...
         real(pyr(bandInds,clr)); zeros(prod(fakePind(size(fakePind,1),:)),1);];
      ch = reconSFpyr(fakePyr, fakePind, [1]);     % recon ori bands only
      [Nly, Nlx] = size(ch);
      im(1:Nly,1:Nlx,clr) = real(expand(im(1:Nly/2,1:Nlx/2,clr),2))/4;
      im(1:Nly,1:Nlx,clr) = im(1:Nly,1:Nlx,clr) + ch;
      vari =  acr0(la+1:la+1,la+1:la+1,nsc,clr);
if cmask(2),
      if vari/D(clr,clr) > 1e-3,
         [im(1:Nly,1:Nlx,clr), snr2(niter,nsc,clr)] = ...
            modacor22(im(1:Nly,1:Nlx,clr), acr0(la-le+1:la+le+1,la-le+1:la+le+1,nsc,clr), p);
      else
        im(1:Nly,1:Nlx,clr) = im(1:Nly,1:Nlx,clr)*sqrt(vari/var2(im(1:Nly,1:Nlx,clr)));
      end
end     % cmask(2)
      im = real(im);

if cmask(1),
      %% Fix marginal stats
      if vari/D(clr,clr) > 1e-3,
          [im(1:Nly,1:Nlx,clr),snr7(niter,2*nsc-1,clr)] = ...
		modskew(im(1:Nly,1:Nlx,clr),skew0p(nsc,clr),p);       % Adjusts skewness
          [im(1:Nly,1:Nlx,clr),snr7(niter,2*nsc,clr)] = ...
		modkurt(im(1:Nly,1:Nlx,clr),kurt0p(nsc,clr),p);         % Adjusts kurtosis
      end
end     % cmask(1)

    end  % for clr = 1:Nclr

  end  %END Coarse-to-fine loop

  for clr = 1:Nclr,

    %% Adjust variance in HP, if higher than desired
if (cmask(2)|cmask(3)|cmask(4)),
    ind = pyrBandIndices(pind,1);
    ch = pyr(ind,clr);
    vHPR = mean2(ch.^2);
    if vHPR > vHPR0(clr),
        ch = ch * sqrt(vHPR0(clr)/vHPR);
        pyr(ind,clr) = ch;
    end
end % cmask(2)


    im(:,:,clr) = im(:,:,clr) + reconSFpyr(real(pyr(:,clr)), pind, [0]);  %recon hi only

if cmask(2),
    la = (Na - 1)/2;
    le = la;
    [im(:,:,clr), snr2(niter,Nsc+2,clr)] = ...
            modacor22(im(:,:,clr), acr0(la-le+1:la+le+1,la-le+1:la+le+1,Nsc+2,clr), p);
end     % cmask(2)
    im = real(im);

    %% Pixel statistics of the PCA channels
    means = mean2(im(:,:,clr));
    vars = var2(im(:,:,clr), means);
if cmask(1),
    im(:,:,clr) = im(:,:,clr) - means;                    % Adjusts mean and variance
    im(:,:,clr) = im(:,:,clr)/sqrt(vars);
    im(:,:,clr) = modskew(im(:,:,clr),skewP(clr)); % Adjusts skewness (keep mean and variance)
    im(:,:,clr) = modkurt(im(:,:,clr),kurtP(clr)); % Adjusts kurtosis (keep mean and variance,
end % if cmask(1)

  end % clr

  %impose desired correlation between RGB bands
  im = reshape(im, Ny*Nx, Nclr);
  im = im - ones(Ny*Nx,Nclr)*diag(mean(im));
  im = adjustCorr1s(im, eye(Nclr));
  im = real(im);
  im = im * sqrt(D) * V.';
  im = im + ones(Ny*Nx,Nclr)*diag(mean0);
  im = reshape(im, Ny, Nx, Nclr);

for clr =1:Nclr,

  %% Pixel statistics of the RGB channels
  means = mean2(im(:,:,clr));
  vars = var2(im(:,:,clr), means);
  [mns mxs] = range2(im(:,:,clr));
  im(:,:,clr) = im(:,:,clr)-means;               % Adjusts mean and variance
  snr7(niter,2*(Nsc+1)+1,clr) = snr(mx0(clr)-mn0(clr),sqrt((mx0(clr)-mxs)^2+(mn0(clr)-mns)^2));
  skews = skew2(im(:,:,clr));
  snr7(niter,2*(Nsc+1)+2,clr) = snr(skew0(clr),skew0(clr)-skews);
  kurts = kurt2(im(:,:,clr));
  snr7(niter,2*(Nsc+1)+3,clr) = snr(kurt0(clr),kurt0(clr)-kurts);
if cmask(1),
  im(:,:,clr) = im(:,:,clr)*sqrt(((1-p)*vars + p*var0(clr))/vars);
end     % cmaks(1)
  im(:,:,clr) = im(:,:,clr)+mean0(clr);
if cmask(1),
  im(:,:,clr) = modskew(im(:,:,clr),skew0(clr)); % Adjusts skewness (keep mean and variance)
  im(:,:,clr) = modkurt(im(:,:,clr),kurt0(clr)); % Adjusts kurtosis (keep mean and variance,
                                                        % but not skewness)
  im(:,:,clr) = max(min(im(:,:,clr),mx0(clr)),mn0(clr)); % Adjusts range (affects everything)
end % if cmask(1)

  %% Force pixels specified by image mask
  if (~isempty(imaskInd))
    im(:,:,clr) = imMask.*imOrig(:,:,clr) + (1-imMask).*im(:,:,clr);
  end
  snr6(niter,clr) = snr(im(:,:,clr)-mean0(clr),im(:,:,clr)-prev_im(:,:,clr));

end  % for clr

  if floor(log2(niter))==log2(niter),
        nq = nq + 1;
        imS(:,:,:,nq) = im;
  end

  tmp = prev_im;
  prev_im=im;
  change = im - tmp;
  
if cmask(5),
  figure(imf);
  subplot(1,2,1);
  showIm(im(:,:,1),[0 255],256/max(Ny,Nx));
  Mach = max(max(max(change)));
  Mich = min(min(min(change)));
  image(uint8(255*(change-Mich)/(Mach-Mich)));axis('image');axis('off');
  title('Change');
  subplot(1,2,2);
  showIm(im(:,:,1),[0 255],256/max(Ny,Nx));
  image(uint8(im)); axis('image');axis('off');
  title(sprintf('iteration %d',niter));
  drawnow
end % cmask(5)

  % accelerator
%  alpha = 0.8;
%  im = im + alpha*change;


end %END  MAIN LOOP

im = prev_im;
im = uint8(im);


function [params] = textureColorAnalysis(im0, Nsc, Nor, Na)

% Analyze texture for application of Portilla-Simoncelli model/algorithm.
%
% [params] = textureColorAnalysis(im0, Nsc, Nor, Na);
%       im0:    original image (uint8, true color)
%       Nsc:    number of scales
%       Nor:    number of orientations
%       Na:     spatial neighborhood considered (Na x Na)
%
% Example: Nsc=4; Nor=4; Na=7;
%
% See also textureColorSynthesis.
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

[Ny,Nx,Nclr] = size(im0);
im0 = double(im0);

% Apply PCA to the RGB components

imPCA = reshape(im0,Ny*Nx,Nclr);
mean0 = mean(imPCA).';
imPCA = imPCA - ones(Ny*Nx,Nclr)*diag(mean0);
Cclr0 = innerProd(imPCA)/(Ny*Nx);
[V,D] = eig(Cclr0);
imPCA = imPCA*V*pinv(sqrt(D));
imPCA = reshape(imPCA,Ny,Nx,Nclr);

%% Check required args are passed
if (nargin < 4)
  error('Function called with too few input arguments');
end

if ( mod(Na,2) == 0 )
  error('Na is not an odd integer');
end

%% If the spatial neighborhood Na is too big for the lower scales,
%% "modacor22.m" will make it as big as the spatial support at
%% each scale:

nth = log2(min(Ny,Nx)/Na);
if nth<Nsc,
  fprintf(1,'Warning: Na will be cut off for levels above #%d !\n', floor(nth+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

la = floor((Na-1)/2);

%% Pixel statistics
mn0 = zeros(Nclr,1);mx0=mn0;var0=mn0;skew0=mn0;kurt0=mn0;
for clr = 1:Nclr,
        [mn0(clr) mx0(clr)] = range2(im0(:,:,clr));
        var0(clr) = var2(im0(:,:,clr));
        skew0(clr) = skew2(im0(:,:,clr));
        kurt0(clr) = kurt2(im0(:,:,clr));
end
statg0 = [mean0 var0 skew0 kurt0 mn0 mx0];

%% "Pixel" statistics of the PCA channels
mnP = zeros(Nclr,1);mxP=mnP;skewP=mnP;kurtP=mnP;
for clr = 1:Nclr,
        [mnP(clr) mxP(clr)] = range2(imPCA(:,:,clr));
        skewP(clr) = skew2(imPCA(:,:,clr));
        kurtP(clr) = kurt2(imPCA(:,:,clr));
end
statgP = [skewP kurtP mnP mxP];

%% Central autoCorr of the PCA bands
acr = NaN * ones(Na,Na,Nsc+2,Nclr);
la = (Na - 1)/2;
for clr = 1:Nclr,
  ac = fftshift(real(ifft2(abs(fft2(imPCA(:,:,clr))).^2)))/(Ny*Nx);
  cy = Ny/2 + 1; cx = Nx/2 + 1;
  le = la;
  ac = ac(cy-le:cy+le,cx-le:cx+le);
  acr(la-le+1:la+le+1,la-le+1:la+le+1,Nsc+2,clr) = ac;
end

%% Build the steerable pyramid
for clr = 1:Nclr,
        [pyr0(:,clr),pind0] = buildSCFpyr(imPCA(:,:,clr),Nsc,Nor-1);
end

if ( any(vectify(mod(pind0,2))) )
  error('Algorithm will fail: Some bands have odd dimensions!');
end

%% Subtract mean of lowBand:
nband = size(pind0,1);
for clr = 1:Nclr,
        pyr0(pyrBandIndices(pind0,nband),clr) = ...
                real(pyr0(pyrBandIndices(pind0,nband),clr))-...
                mean(real(pyr0(pyrBandIndices(pind0,nband),clr)));
end
rpyr0 = real(pyr0);
apyr0 = abs(pyr0);

figure(gcf)
showIm(im0(:,:,1),[0 255],1);
image(uint8(im0)); axis('image');axis('off');title('Original');  drawnow

%% Subtract mean of magnitude:
magMeans0 = zeros(size(pind0,1), Nclr);
for clr = 1:Nclr,
for nband = 1:size(pind0,1)
  indices = pyrBandIndices(pind0,nband);
  magMeans0(nband,clr) = mean(apyr0(indices,clr));
  apyr0(indices,clr) = apyr0(indices,clr) - magMeans0(nband,clr);
end
end

%% Compute central autoCorr of lowband
im = zeros(Ny,Nx,Nclr);
skew0p = zeros(Nsc+1,Nclr);
kurt0p = zeros(Nsc+1,Nclr);
for clr = 1:Nclr,
  nband = size(pind0,1);
  ch = pyrBand(pyr0(:,clr),pind0,nband);
  [Nly Nlx] = size(ch);
  [mpyr,mpind] = buildSFpyr(real(ch),0,0);
  im(1:Nly,1:Nlx,clr) = pyrBand(mpyr,mpind,2);
  Sch = min(Nly,Nlx); %size of low band
  le = min(Sch/2-1,la);
  cy = Nly/2+1;
  cx = Nlx/2+1;
  ac = fftshift(real(ifft2(abs(fft2(im(1:Nly,1:Nlx,clr))).^2)))/prod(size(ch));
  ac = ac(cy-le:cy+le,cx-le:cx+le);
  acr(la-le+1:la+le+1,la-le+1:la+le+1,Nsc+1,clr) = ac;
  vari = ac(le+1,le+1);
  if vari/D(clr,clr) > 1e-6,
        skew0p(Nsc+1,clr) = mean2(im(1:Nly,1:Nlx,clr).^3)/vari^1.5;
        kurt0p(Nsc+1,clr) = mean2(im(1:Nly,1:Nlx,clr).^4)/vari^2;
  else
        skew0p(Nsc+1,clr) = 0;
        kurt0p(Nsc+1,clr) = 3;
  end
end

%% Compute  central autoCorr of each Mag band, and the autoCorr of the
%% combined (non-oriented) band.
ace = NaN * ones(Na,Na,Nsc,Nor,Nclr);
for clr = 1:Nclr,
for nsc = Nsc:-1:1,
  for nor = 1:Nor,
    nband = (nsc-1)*Nor+nor+1;
    ch = pyrBand(apyr0(:,clr),pind0,nband);
    [Nly, Nlx] = size(ch);
    Sch = min(Nlx, Nly);
    le = min(Sch/2-1,la);
    cx = Nlx/2+1;  %Assumes Nlx even
    cy = Nly/2+1;
    ac = fftshift(real(ifft2(abs(fft2(ch)).^2)))/prod(size(ch));
    ac = ac(cy-le:cy+le,cx-le:cx+le);
    ace(la-le+1:la+le+1,la-le+1:la+le+1,nsc,nor,clr) = ac;
  end

  %% Combine ori bands

  bandNums = [1:Nor] + (nsc-1)*Nor+1;  %ori bands only
  ind1 = pyrBandIndices(pind0, bandNums(1));
  indN = pyrBandIndices(pind0, bandNums(Nor));
  bandInds = [ind1(1):indN(length(indN))];
  %% Make fake pyramid, containing dummy hi, ori, lo
  fakePind = [pind0(bandNums(1),:);pind0(bandNums(1):bandNums(Nor)+1,:)];
  fakePyr = [zeros(prod(fakePind(1,:)),1);...
         rpyr0(bandInds,clr); zeros(prod(fakePind(size(fakePind,1),:)),1);];
  ch = reconSFpyr(fakePyr, fakePind, [1]);     % recon ori bands only
  im(1:Nly,1:Nlx,clr) = real(expand(im(1:Nly/2,1:Nlx/2,clr),2))/4;
  im(1:Nly,1:Nlx,clr) = im(1:Nly,1:Nlx,clr) + ch;
  ac = fftshift(real(ifft2(abs(fft2(im(1:Nly,1:Nlx,clr))).^2)))/prod(size(ch));
  ac = ac(cy-le:cy+le,cx-le:cx+le);
  acr(la-le+1:la+le+1,la-le+1:la+le+1,nsc,clr) = ac;
  vari = ac(le+1,le+1);
  if vari/var0(clr) > 1e-6,
        skew0p(nsc,clr) = mean2(im(1:Nly,1:Nlx,clr).^3)/vari^1.5;
        kurt0p(nsc,clr) = mean2(im(1:Nly,1:Nlx,clr).^4)/vari^2;
  else
        skew0p(nsc,clr) = 0;
        kurt0p(nsc,clr) = 3;
  end
end
end

%% Compute the cross-correlation matrices of the coefficient magnitudes
%% pyramid at the different levels and orientations

C0 = zeros(Nclr*Nor,Nclr*Nor);
Cx0 = zeros(Nclr*Nor,Nclr*Nor,Nsc-1);
Cr0 = zeros(2*Nor*Nclr,Nclr*2*Nor,Nsc+1);
Crx0 = zeros(Nclr*Nor,2*Nclr*Nor,Nsc);

for nsc = 1:Nsc,
  firstBnum = (nsc-1)*Nor+2;
  cousinSz = prod(pind0(firstBnum,:));
  ind = pyrBandIndices(pind0,firstBnum);
  cousinInd = ind(1) + [0:Nor*cousinSz-1];

  cousins = zeros(cousinSz,Nor,Nclr);
  rcousins = zeros(cousinSz,Nor,Nclr);
  parents = zeros(cousinSz,Nor,Nclr);
  rparents = zeros(cousinSz,2*Nor,Nclr);

  for clr = 1:Nclr,

  if (nsc<Nsc)
    for nor=1:Nor,
      nband = (nsc-1+1)*Nor+nor+1;

      tmp = expand(pyrBand(pyr0(:,clr), pind0, nband),2)/4;
      rtmp = real(tmp); itmp = imag(tmp);
      %% Double phase:
      tmp = sqrt(rtmp.^2 + itmp.^2) .* exp(2 * sqrt(-1) * atan2(rtmp,itmp));
      rparents(:,nor,clr) = vectify(real(tmp));
      rparents(:,Nor+nor,clr) = vectify(imag(tmp));

      tmp = abs(tmp);
      parents(:,nor,clr) = vectify(tmp - mean2(tmp));
    end
  else
    tmp = real(expand(pyrLow(pyr0(:,clr),pind0),2))/4;
    rparents(:,1:5,clr) = [vectify(tmp),...
                vectify(shift(tmp,[0 2])), vectify(shift(tmp,[0 -2])), ...
                vectify(shift(tmp,[2 0])), vectify(shift(tmp,[-2 0]))];

    parents = [];
  end

  cousins(:,:,clr) = reshape(apyr0(cousinInd,clr), [cousinSz Nor]);
  rcousins(:,:,clr) = reshape(real(pyr0(cousinInd,clr)), [cousinSz Nor]);

  end % clr

  nc = size(cousins,2)*Nclr;   np = size(parents,2)*Nclr;
  cousins = reshape(cousins,[cousinSz nc]);
  parents = reshape(parents,[cousinSz np]);
  C0(:,:,nsc) = innerProd(cousins)/cousinSz;
  if (np > 0)
    Cx0(:,:,nsc) = (cousins'*parents)/cousinSz;
  end

  if (nsc == Nsc),
     rparents = rparents(:,1:5,:);
  end   
  nrp = size(rparents,2)*Nclr;
  nrc = size(rcousins,2)*Nclr;   
  rcousins = reshape(rcousins,[cousinSz nrc]);
  rparents = reshape(rparents,[cousinSz nrp]);

  Cr0(1:nrc,1:nrc,nsc) = innerProd(rcousins)/cousinSz;
  if (nrp > 0)
    Crx0(1:nrc,1:nrp,nsc) = (rcousins'*rparents)/cousinSz;
    if (nsc==Nsc)
      Cr0(1:nrp,1:nrp,Nsc+1) = innerProd(rparents)/cousinSz;
    end
  end
end


%% Calculate the variance of the HF residual.

vHPR0 = zeros(Nclr,1);
for clr = 1:Nclr,
        channel = pyr0(pyrBandIndices(pind0,1),clr);
        vHPR0(clr) = mean2(channel.^2);
end

statsLPim = [skew0p; kurt0p];

params = struct('pixelStats', statg0, ...
                'pixelStatsPCA', statgP, ...
        		'pixelLPStats', statsLPim, ...
                'autoCorrReal', acr, ...
                'autoCorrMag', ace, ...
                'magMeans', magMeans0, ...
                'cousinMagCorr', C0, ...
                'parentMagCorr', Cx0, ...
                'cousinRealCorr', Cr0, ...
                'parentRealCorr', Crx0, ...
                'varianceHPR', vHPR0,...
                'colorCorr', Cclr0);



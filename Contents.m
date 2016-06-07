% Tools for analysis and synthesis of texture images.
% Version 1.0 March 2001.
%
% Created: Javier Portilla and Eero Simoncelli
%          javier@decsai.ugr.es / eero.simoncelli@nyu.edu
%
% See Readme.txt file for a brief description.
% See ChangeLog file for latest modifications. 
% Type "help <command-name>" for documentation on individual commands.
% -----------------------------------------------------------------
% Demonstrations: 
%   example1	     - synthesis of "random text"
%   example2         - synthesis mixing real and synthetic data
%
% Primary entry points:
%   textureAnalysis  - Extract a set of parameters for a texture image
%   textureSynthesis - Synthesize a new texture from a set of parameters.
%
% Utility functions:
%   snr       - compute signal-to-noise ratio in dB
%   expand    - resample at higher resolution (in Fourier domain)
%   shrink    - resample at lower resolution (in Fourier domain)
%
% Functions that perform projection parameter constraint surfaces:
%   modacor22	 - modify autocorrelation of an image
%   modkurt	 - modify kurtosis (4th moment divided by squared variance)
%   modskew	 - modify skewness (3rd moment divided by variance^1.5)
%   adjustCorr1s - modify cross-correlation
%   adjustCorr2s - modify cross-correlation, with some variables held fixed.
%	
% Example texture images:	
%   text.pgm			(scanned by us)
%   nuts.pgm			(from VisTex texture data base)
%   metal.pgm			(VisTex)
%   reptil_skin.pgm             (Brodatz)
%   checkerboard.pgm		(artificial)
%   sawtooth.pgm		(artificial)

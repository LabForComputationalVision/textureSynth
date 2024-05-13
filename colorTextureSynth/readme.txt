=====================================================
Color texture analysis & synthesis toolbox for Matlab

J. Portilla and E. P. Simoncelli
portilla@io.cfmac.csic.es, eero.simoncelli@nyu.edu

October 2000, New York University, New York
Released Version 1.0: January 2009, CSIC, Madrid.

=====================================================

This toolbox implements a color version of the texture model and
synthesis method described in:

   J Portilla and E P Simoncelli, "A parametric texture model
   based on joint statistics of complex wavelet coefficients",
   International Journal of Computer Vision, 40(1):49-71, 2000

Instructions:

1) You'll need to first download the original (grayscale) toolbox from
http://www.cns.nyu.edu/~eero/texture/

2) Save to a folder and make that folder available in your Matlab
path, eg.  addpath(genpath(<myFolder>)) - see help information for
"path" in Matlab.

3) add the files "textureColorAnalysis.m" and
"textureColorSynthesis.m" to that folder

4) run the file "demo.m" script, which will generate a synthetic
texture based on the example image provided (olives256.o.bmp).
The script is easily modified to run on any image file you like.

Please see copyright statement in the "copyright.txt" file and report
any bugs to  "portilla@io.cfmac.csic.es".

Enjoy!

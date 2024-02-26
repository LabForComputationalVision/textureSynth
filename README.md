# textureSynth

This package contains MatLab code for analyzing and synthesizing digital
image of visual texture. The algorithm is described in the references
given at the bottom of this document. Further information, as well the
most recent versions of the code, are available at

    http://www.cns.nyu.edu/~lcv/texture/

Incremental changes to the code are documented in the `ChangeLog` file.

Written by Javier Portilla and Eero Simoncelli, 1999-2000.

Comments/Suggestions/Bugs to:
- javier\@decsai.ugr.es 
- eero.simoncelli\@nyu.edu

## INSTALLATION

1)  download and unpack the code. You can put the code anywhere on your
    system, but we'll assume it's in a directory named `textureSynth`.

2)  download and unpack the matlabPyrTools package from
    https://github.com/LabForComputationalVision/matlabPyrTools This is
    a collection of tools for multi-scale decomposition of images. You
    can put the code anywhere on your system, but we'll assume it's in a
    directory named `matlabPyrTools`. Please use version 1.4 or newer of
    the matlabPyrTools.

3)  Run matlab, and put the matlabPyrTools directory in your path: >
    path('matlabPyrTools', path);

4)  The matlabPyrTools distribution includes a MEX subdirectory
    containing binary executables, precompiled for various platforms
    (SunOS,Solaris, Linux,Windows). You may need to recompile these on
    your platform. In addition, you should either move the relavent
    files from the MEX subdirectory into the main directory, OR create a
    link/alias to them, OR place the MEX subdirectory in your matlab
    path.

## USING THE SOFTWARE

-   To see a demonstration, start Matlab, change directories (using
    `cd`) into the textureSynth directory, and execute `example1`.

-   If you want to learn how to use the texture analysis and synthesis
    functions, take a look at example1.m and example2.m (these include
    many explanatory comments).

-   For a listing of matlab function included in this package, execute `help
    textureSynth`

-   For details on any of the functions, execute `help <name_of_function>`

## REFERENCES

J Portilla and E P Simoncelli. A Parametric Texture Model based on Joint
Statistics of Complex Wavelet Coefficients. Int'l Journal of Computer
Vision. 40(1):49-71, October, 2000.
http://www.cns.nyu.edu/\~eero/ABSTRACTS/portilla99-abstract.html

J Portilla and E P Simoncelli Texture Modeling and Synthesis using Joint
Statistics of Complex Wavelet Coefficients. IEEE Workshop on Statistical
and Computational Theories of Vision, Fort Collins, CO, 22 June 1999.\
http://www.cns.nyu.edu/\~eero/ABSTRACTS/portilla99a-abstract.html

J Portilla and E P Simoncelli. Texture Representation and Synthesis
Using Correlation of Complex Wavelet Coefficient Magnitudes. Technical
Report #54, Consejo Superior de Investigaciones Cientificas (CSIC),
Madrid. 29 March 1999.\
http://www.cns.nyu.edu/\~eero/ABSTRACTS/portilla98-abstract.html

E P Simoncelli and J Portilla. Texture Characterization via Joint
Statistics of Wavelet Coefficient Magnitudes. In 5th IEEE Int'l Conf on
Image Processing. Chicago, IL. Oct 4-7, 1998.\
http://www.cns.nyu.edu/\~eero/ABSTRACTS/simoncelli98b-abstract.html

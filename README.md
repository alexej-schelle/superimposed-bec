# superimposed-quantum-fields

The routine superimposed-quantum-fields.py is a Markov sampling method, more specific, a (generalised) Metropolis sampling algorithm for calculating the integrated (random) superimposed condensate and non-condensate wave fields (order parameters) for a given total particle number, temperature and trap frequency of a Bose-Einstein condensate.

For proper installation, please replace the path /your-installation-path/ with your installation path in the lines ... of this code.

This software version corresponds to the publication Fluctuation and Noise Letters, Vol. 16, No. 01, 1750009 (2017).

The Python routine superimposed-quantum-fields.py is compatible with Python 2 or Python 3 (for Version 3, some parts have to be changed) 
and calculates :

- Condensate field modes
- Condensate wave fields at equilibrium (symmetry aspects)
- Condensate wave field propagations (from low energy to high energy)
- Partial phase distributions of real valued field modes
- Partial phase distributions of imaginary valued field modes
- Total phases of the condensate wave field
- Chemical potentials of the condensate

In order to start bec_symmytry_breaking.py, please specify BEC parameters :

maxmode # typical mode size for analysis : 50 - 500 modes #
ptn # typical particle number : 10^3 - 10^5 #
sample # Typical sample size : 10^5 - 10^8 #

omx # trap frequency in x direction #
omy # trap frequency in y direction #
omz # trap frequency in z direction #

start_temp # in units of nK #

in the source code as required and run the simulation in the command line with the command

python bec_symmytry_breaking.py

using Python Version 2 or Version 3 (for Version 3, some parts in the code must be changed).

The source code is granted with MIT LICENSE (Copyright, 2017) :

-1- :

License Copyright : Dr. A. Schelle, Bachschmidstr. 4, 87600 Kaufbeuren
License Type : MIT License (2022)
License Contact : E-Mail : alexej.schelle@gmail.com

-2- :

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the Software "bec_symmetry_breaking.py"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice -1- and this permission notice -2- shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

-3- : Free support is provided per email at support@krealix.de.

-4- : For more explicit consulting and discussions, please schedule a meeting at https://calendly.com/alexej-schelle/.

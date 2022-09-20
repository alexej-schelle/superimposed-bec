# superimposed-quantum-fields

The routine superimposed-bec.py is a Markov sampling method, more specific, a (generalised) Metropolis sampling algorithm for calculating the integrated (random) superimposed condensate and non-condensate wave fields (order parameters) for a given total particle number, temperature and trap frequency of a Bose-Einstein condensate, as it would correspond to a measurement of the two orthogonal quantum fields (condensate and non-condensate).

For proper installation, please replace the path /your-installation-path/ with your installation path in the lines ... of this code.

The Python routine superimposed-quantum-fields.py is compatible with Python 2 or Python 3 (for Version 3, some parts have to be changed) 
and calculates :

- Superimposed field modes

In order to start superimposed-quantum-fields.py, please specify BEC parameters :

maxmode # typical mode size for analysis : 50 - 500 modes # <br>
ptn # typical particle number : 10^3 - 10^5 # <br>
sample # Typical sample size : 10^5 - 10^8 # <br>

omx # trap frequency in x direction # <br>
omy # trap frequency in y direction # <br>
omz # trap frequency in z direction # <br>

start_temp # in units of nK # <br>

in the source code as required and run the simulation in the command line with the command

python superimposed-quantum-fields.py

using Python Version 2 or Version 3 (for Version 3, some parts in the code must be changed).

The source code is granted with MIT LICENSE (Copyright, 2022) :

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

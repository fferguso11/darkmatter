# darkmatter
Dark Matter Project

Included is a simulation that produces a highly-populated recoil spectrum for the LUX experiment (the target) for a given WIMP mass and WIMP-nucleon cross section, and then generates sample datasets for a single run of LUX at each point in a 100x100 mass/cross section phase space. A K.S test is run for the target and each sample, with 50 iterations at each point to produce a mean p-value. These p-values help to quantify the ability of these measurements to identify the mass and cross section.

The velocity distribution on the samples can be varied to explore the systematic error introduced by uncertainty in this. 

'yield.dat' and 'formfactor.dat' represent input files whose raw data is used in the main simulation
- will need to be placed in same folder as the simulation when run

Consult the accompanying report for more information

Before running simulation, open the code to define the mass and cross section for the target and specify the output file. The default setting is for no file output, and target Mw = 60, sigma = 4D-44.


LICENCE AGREEMENT:

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation files
(the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

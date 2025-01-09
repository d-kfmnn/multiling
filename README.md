[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# MultiLinG - AIG Verification using Linear Extractions from Groebner Bases


Our tool MultiLinG able to verify AIGs using linear extractions from Groebner bases.

Dependencies: `libgmp` (https://gmplib.org/)
              `msolve` (https://msolve.lip6.fr/)
              

Use 
    `./configure.sh && make` 
    
to configure and build `MultiLinG`.


    usage : multiling <input.aig> <mode | spec_file> [<option> ...] 

    Either a <mode> with an automatically generated specification has to 
    be selected or a user-defined spec has to be provided in a file <spec_file>


    <mode> = -mult:
        assumes that the AIG represents an unsigned n*n = 2n multiplier
        automatically generated spec: 2^(2n-1)s_(2n-1) + ...+s_0 - 
        (2^(n-1)a_(n-1)+...+a_0)*(2^(n-1)b_(n-1)+...+b_0) 

    <mode> = -miter:
        assumes that the AIG represents an unsigned n*n = 2n multiplier
        automatically generated spec: s-1 

    <spec_file> =    file containing a single polynomial used as the specification


    <option> = the following options are available 
      -h | --help           print this command line summary 
      -v<1,2,3,4>           different levels of verbosity(default -v1) 
      -no-counter-examples  do not generate and write counter examples


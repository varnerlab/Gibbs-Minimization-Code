### Network Free Gibbs Minimization
Code to solve the constrained Gibbs energy minimization problem for a system with multiple reactions encoded in the [Julia](https://julialang.org) programming language.

### Installation and Requirements
You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/Gibbs-Minimization-Code

or

	$ git clone https://github.com/varnerlab/Gibbs-Minimization-Code

To execute a Gibbs minimization job, Julia must be installed on your machine along with the [JSON](https://github.com/JuliaIO/JSON.jl) and
[PyPlot](https://github.com/JuliaPy/PyPlot.jl) Julia packages. Julia can be downloaded/installed on any platform.
The required [Julia](http://julialang.org) packages can be installed by executing the commands:

	julia> Pkg.add("PyPlot")

and

	julia> Pkg.add("JSON")

in the Julia REPL.

### Where is the thermodynamic data stored?
The Gibbs energy minimization calculation uses thermodynamic data for different chemical species that is stored in the [Database.json](https://github.com/varnerlab/Gibbs-Minimization-Code/blob/master/data/Database.json) file.
This JSON file has species records of the form:

    {
      "delta_gibbs_in_kj_mol":"-426.71",
      "symbol":"glucose",
      "element_array":{
        "C":"6",
        "H":"12",
        "O":"6",
        "N":"0",
        "P":"0"
      }
    }

which holds the Gibbs energy of formation in the ``delta_gibbs_in_kj_mol`` field, and the species symbol in the ``symbol`` field, and the
molecular composition in the ``element_array`` array.

### How do we run a calculation?
The Gibbs minimization can be run by executing the ``global_gibbs_minimization.jl`` script. In this script you specify the chemical species in your system
(must match the ``symbol`` records in the ``Database.json`` file), the initial composition of the system (in units of mmol/L),
the system volume (in units of L) and the lower and upper bounds for each system species (in units of mmol/L).

The script solves the minimization problem ``number_of_runs`` times, with the best results of the ith calculation serving as the starting point for
calculation i+1. At the end of ``number_of_runs`` calculations, the best solutions (and associated error) are written to the ``best_soln_archive.dat`` and
``best_error_archive.dat`` text files. The species abundance (in units of mmol/L) is written in the first ``number_of_species`` rows
while the Lagrangian multipliers are written in the remaining entries. Each column is a separate potential solution.

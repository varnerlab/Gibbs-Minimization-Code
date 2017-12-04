include("Includes.jl")

# hard constants -
VOLUME = 1.186e-13      # volume of cell L
FACTOR = 1e-12          # convert mmol to fmol


function check_element_balances(soln_array,problem_data_model)

    # get some stuff from the data model -
    atom_matrix = problem_data_model.atom_matrix
    total_atom_array = problem_data_model.total_atom_array
    number_of_species = problem_data_model.number_of_species

    # compute the elemental residual -
    element_error = atom_matrix*soln_array[1:number_of_species] - total_atom_array

    # 5 x 1 array
    return transpose(element_error)*element_error
end

function calculate_mole_fraction(state_array::Array{Float64,1})

    # compute the total moles -
    total_mol = sum(state_array)

    # scale -
    mol_fraction_array = (1.0/total_mol)*state_array

    # return -
    return (mol_fraction_array,total_mol)
end

function objective_function(parameter_guess,problem_data_model)

    # hard constants -
    # VOLUME = 1.186e-13
    # FACTOR = 1e-12

    # get data from the problem object J/fmol
    R_constant = 8.314*(FACTOR/1000)

    # get stuff from the problem object -
    gibbs_energy_array_in_j_mol = problem_data_model.gibbs_energy_array_in_j_mol
    atom_matrix = problem_data_model.atom_matrix
    total_atom_array = problem_data_model.total_atom_array

    # split the parameter array into state and lambda -
    number_of_species = problem_data_model.number_of_species
    state_array = parameter_guess[1:number_of_species]
    lambda_array = parameter_guess[(number_of_species+1):end]

    # what is the system temp?
    system_temperature_in_kelvin = problem_data_model.system_temperature_in_kelvin

    # compute the scaling factor -
    scaling_factor = (R_constant*system_temperature_in_kelvin)

    # compute the mole fraction -
    (mole_fraction_array,total_mol) = calculate_mole_fraction(state_array)

    # compute the energy balances -
    energy_balances = (1/scaling_factor)*gibbs_energy_array_in_j_mol+log.(e,mole_fraction_array)+(1/scaling_factor)*transpose(atom_matrix)*lambda_array

    # compute the mass balances -
    mass_balances = atom_matrix*state_array - total_atom_array

    # compute scale factors -
    norm_energy_balances = norm(energy_balances)
    norm_mass_balances = norm(mass_balances)

    #@show (norm_energy_balances,norm_mass_balances)

    # build error array -
    error_array = [
        energy_balances ; 1000*mass_balances
    ]

    # return -
    return transpose(error_array)*error_array
end



function main(species_list,species_dictionary,
    initial_composition_dictionary,
    lower_bound_array,
    upper_bound_array,
    initial_multiplier_dictionary,
    system_temperature_in_kelvin)

    # hard constants -
    # VOLUME = 1.186e-13
    # FACTOR = 1e-12

    # Build the atom_matrix -
    atom_matrix = buildAtomMatrix("./data/Database.json",species_list)

    # build the energy of formation array -
    gibbs_energy_array_in_j_mol = buildGibbsEnergyOfFormationArray("./data/Database.json",species_list)

    # build a data model -
    problem_data_model = CREProblemDataModel()
    problem_data_model.species_dictionary = species_dictionary
    problem_data_model.initial_composition_dictionary = initial_composition_dictionary
    problem_data_model.species_list = species_list
    problem_data_model.system_temperature_in_kelvin = system_temperature_in_kelvin
    problem_data_model.atom_matrix = atom_matrix
    problem_data_model.gibbs_energy_array_in_j_mol = gibbs_energy_array_in_j_mol
    problem_data_model.number_of_species = length(species_list)

    # calculate the initial atom value -
    tmp_mol_array::Array{Float64,1} = Float64[]
    for (index,species_key) in enumerate(species_list)

        value = initial_composition_dictionary[species_key]
        push!(tmp_mol_array,value)
    end

    # initial atom vector (A)
    A = atom_matrix*tmp_mol_array
    problem_data_model.total_atom_array = A

    # setup the initial condition -
    initial_parameter_guess::Array{Float64,1} = Float64[]
    for key in species_list

        value = initial_composition_dictionary[key]
        push!(initial_parameter_guess,value)
    end

    # C,H,O,N -
    push!(initial_parameter_guess,initial_multiplier_dictionary["C"])
    push!(initial_parameter_guess,initial_multiplier_dictionary["H"])
    push!(initial_parameter_guess,initial_multiplier_dictionary["O"])
    push!(initial_parameter_guess,initial_multiplier_dictionary["N"])
    push!(initial_parameter_guess,initial_multiplier_dictionary["P"])

    # setup the lower and upper bounds -
    problem_data_model.parameter_lower_bound_array = (lower_bound_array)*(VOLUME/FACTOR)
    problem_data_model.parameter_upper_bound_array = (upper_bound_array)*(VOLUME/FACTOR)

    # call the minimizer -
    (error_archive,parameter_archive) = sa_estimate(objective_function,initial_parameter_guess,problem_data_model)

    # check the element error -
    min_index = indmin(error_archive)
    best = parameter_archive[:,min_index]
    element_error = check_element_balances(best,problem_data_model)

    # return -
    return (error_archive,element_error,parameter_archive)
end

# how many atoms are we balancing over?
number_of_elements = 5

# specify the metabolites I want to optimize -
species_list = [
    "glucose"                                               ;   # 1
    "fructose-6-phosphate"                                  ;   # 2
    "glucose-6-phosphate"                                   ;   # 3
    "atp"                                                   ;   # 4
    "adp"                                                   ;   # 5
    "6-phospho-D-gluconate"                                 ;   # 6
    "D-ribulose-5-phosphate"                                ;   # 7
    "5-phospho-alpha-D-ribose-1-diphosphate"                ;   # 8
    "D-xylulose-5-phosphate"                                ;   # 9
    "nicotinamide-adenine-dinucleotide-phosphate"           ;   # 10
    "nicotinamide-adenine-dinucleotide-phosphate-reduced"   ;   # 11
    "co2"                                                   ;   # 12
    "D-xylulose-5-phosphate"                                ;   # 13
    "alpha-D-Ribose-5-phosphate"                            ;   # 14
    "sedoheptulose-7-phosphate"                             ;   # 15
    "D-erythrose-4-phosphate"                               ;   # 16
    "D-fructose-1,6-bisphosphate"                           ;   # 17
    "dihydroxyacetone-phosphate"                            ;   # 18
    "3-phospho-D-glyceroyl-phosphate"                       ;   # 19
    "3-phospho-D-glycerate"                                 ;   # 20
    "D-glycerate-2-phosphate"                               ;   # 21
    "phosphoenolpyruvate"                                   ;   # 22
    "pyruvate"                                              ;   # 23
    "D-lactate"                                             ;   # 24
    "inorganic-phosphate"                                   ;   # 25
    "amp"                                                   ;   # 26
    "water"                                                 ;   # 27
    "glyceraldehyde-3-phosphate"                            ;   # 28
]

# number of states?
number_of_species = length(species_list)

# Build the species dictionary -
species_dictionary = buildSpeciesDictionary("./data/Database.json")
atom_matrix = buildAtomMatrix("./data/Database.json",species_list)

# set the T -
system_temperature_in_kelvin = 298.15

# ok, setup a number of hot starts -
best_error_archive = Float64[]
element_error_archive = Float64[]
best_soln_archive = zeros(number_of_species+number_of_elements,1)
number_of_runs = 100
for run_index = 1:number_of_runs

    # initialize -
    initial_composition_dictionary = Dict()
    initial_multiplier_dictionary = Dict()

    # do we have a temp.out file -
    if (isfile("temp.out") == true)

        # open the file -
        tmp_data = readdlm("temp.out")

        # grab the state data -
        state_array = tmp_data[1:number_of_species]
        for (index,species_symbol) in enumerate(species_list)
            initial_composition_dictionary[species_symbol] = state_array[index]
        end

        # grab the multiplier data -
        lambda_array = tmp_data[(number_of_species+1):end]
        initial_multiplier_dictionary["C"] = lambda_array[1]
        initial_multiplier_dictionary["H"] = lambda_array[2]
        initial_multiplier_dictionary["O"] = lambda_array[3]
        initial_multiplier_dictionary["N"] = lambda_array[4]
        initial_multiplier_dictionary["P"] = lambda_array[5]
    else

        # setup the composition -
        initial_composition_dictionary["glucose"] = 1.0*(VOLUME/FACTOR)
        initial_composition_dictionary["atp"] = 0.034*(VOLUME/FACTOR)
        initial_composition_dictionary["adp"] = 6.4*(VOLUME/FACTOR)
        initial_composition_dictionary["glucose-6-phosphate"] =  0.041*(VOLUME/FACTOR)
        initial_composition_dictionary["fructose-6-phosphate"] = 0.016*(VOLUME/FACTOR)
        initial_composition_dictionary["6-phospho-D-gluconate"] =  0.026*(VOLUME/FACTOR)
        initial_composition_dictionary["D-ribulose-5-phosphate"] = 0.087*(VOLUME/FACTOR)
        initial_composition_dictionary["5-phospho-alpha-D-ribose-1-diphosphate"] = 0.026*(VOLUME/FACTOR)
        initial_composition_dictionary["D-xylulose-5-phosphate"] =  0.013*(VOLUME/FACTOR)
        initial_composition_dictionary["nicotinamide-adenine-dinucleotide-phosphate"] = 0.075*(VOLUME/FACTOR)
        initial_composition_dictionary["nicotinamide-adenine-dinucleotide-phosphate-reduced"] = 0.075*(VOLUME/FACTOR)
        initial_composition_dictionary["co2"] = 1.2*(VOLUME/FACTOR)
        initial_composition_dictionary["alpha-D-Ribose-5-phosphate"] = 0.078*(VOLUME/FACTOR)
        initial_composition_dictionary["sedoheptulose-7-phosphate"] =  0.016*(VOLUME/FACTOR)
        initial_composition_dictionary["D-erythrose-4-phosphate"] = 0.064*(VOLUME/FACTOR)
        initial_composition_dictionary["D-fructose-1,6-bisphosphate"] = 0.146*(VOLUME/FACTOR)
        initial_composition_dictionary["dihydroxyacetone-phosphate"] = 0.01*(VOLUME/FACTOR)
        initial_composition_dictionary["3-phospho-D-glyceroyl-phosphate"] = 0.023*(VOLUME/FACTOR)
        initial_composition_dictionary["3-phospho-D-glycerate"] =  0.01*(VOLUME/FACTOR)
        initial_composition_dictionary["D-glycerate-2-phosphate"] = 0.084*(VOLUME/FACTOR)
        initial_composition_dictionary["phosphoenolpyruvate"] = 0.011*(VOLUME/FACTOR)
        initial_composition_dictionary["pyruvate"] = 0.234*(VOLUME/FACTOR)
        initial_composition_dictionary["inorganic-phosphate"] = 1.00*(VOLUME/FACTOR)
        initial_composition_dictionary["D-lactate"] =  1.7*(VOLUME/FACTOR)
        initial_composition_dictionary["amp"] = 0.01*(VOLUME/FACTOR)
        initial_composition_dictionary["water"] =  0.01*(VOLUME/FACTOR)
        initial_composition_dictionary["D-xylulose-5-phosphate"] = 0.01*(VOLUME/FACTOR)
        initial_composition_dictionary["glyceraldehyde-3-phosphate"] = 0.00405*(VOLUME/FACTOR)

        # setup multiplier -
        initial_multiplier_dictionary["C"] = 100.0
        initial_multiplier_dictionary["H"] = 100.0
        initial_multiplier_dictionary["O"] = 100.0
        initial_multiplier_dictionary["N"] = 100.0
        initial_multiplier_dictionary["P"] = 100.0
    end


    # Setup the bounds -
    lower_bound_array = [

        0.0001  ;    # "glucose"                                               ;   # 1
        0.00001 ;   # "fructose-6-phosphate"                                  ;   # 2
        0.00001 ;     # "glucose-6-phosphate"                                   ;   # 3
        0.00001 ;     # "atp"                                                   ;   # 4
        0.00001 ;     # "adp"                                                   ;   # 5

        0.00001 ;     # "6-phospho-D-gluconate"                                 ;   # 6
        0.00001 ;     # "D-ribulose 5-phosphate"                                ;   # 7
        0.00001 ;     # "5-phospho-alpha-D-ribose-1-diphosphate"                ;   # 8
        0.00001 ;     # "D-xylulose-5-phosphate"                                ;   # 9
        1.0     ;    # "nicotinamide-adenine-dinucleotide-phosphate"           ;   # 10

        1.2     ;    # "nicotinamide-adenine-dinucleotide-phosphate-reduced"   ;   # 11
        0.0001  ;     # "co2"                                                   ;   # 12
        0.00001 ;     # "D-xylulose-5-phosphate"                                ;   # 13
        0.00001 ;     # "alpha-D-Ribose-5-phosphate"                            ;   # 14
        0.00001 ;     # "sedoheptulose-7-phosphate"                             ;   # 15

        0.00001 ;     # "D-erythrose-4-phosphate"                               ;   # 16
        0.00001 ;     # "D-fructose-1,6-bisphosphate"                           ;   # 17
        0.00001 ;     # "dihydroxyacetone-phosphate"                            ;   # 18
        0.00001 ;     # "3-phospho-D-glyceroyl-phosphate"                       ;   # 19
        0.00001 ;     # "3-phospho-D-glycerate"                                 ;   # 20

        0.00001 ;     # "D-glycerate-2-phosphate"                               ;   # 21
        0.00001 ;     # "phosphoenolpyruvate"                                   ;   # 22
        0.00001 ;     # "pyruvate"                                              ;   # 23
        0.00001 ;     # "D-lactate"                                             ;   # 24
        0.0001 ;     # "inorganic-phosphate"                                   ;   # 25
        0.00001 ;     # "amp"                                                   ;   # 26
        0.00001 ;     # "water"                                                 ;   # 27
        0.00001 ;     # "glyceraldehyde-3-phosphate"                            ;   # 28

        -Inf ;    # l-C
        -Inf ;    # l-H
        -Inf ;    # l-O
        -Inf ;    # l-N
        -Inf ;    # l-P
    ]

    upper_bound_array = [

        100.0 ;       # "glucose"                                               ;   # 1
        10.0 ;       # "fructose-6-phosphate"                                  ;   # 2
        10.0 ;       # "glucose-6-phosphate"                                   ;   # 3
        10.0 ;       # "atp"                                                   ;   # 4
        10.0 ;       # "adp"                                                   ;   # 5
        10.0 ;       # "6-phospho-D-gluconate"                                 ;   # 6
        10.0 ;       # "D-ribulose-5-phosphate"                                ;   # 7
        10.0 ;       # "5-phospho-alpha-D-ribose-1-diphosphate"                ;   # 8
        10.0 ;       # "D-xylulose-5-phosphate"                                ;   # 9
        150.0 ;     # "nicotinamide-adenine-dinucleotide-phosphate"           ;   # 10

        150.0 ;     # "nicotinamide-adenine-dinucleotide-phosphate-reduced"   ;   # 11
        25.0 ;       # "co2"                                                   ;   # 12
        10.0 ;       # "D-xylulose-5-phosphate"                                ;   # 13
        10.0 ;       # "alpha-D-Ribose-5-phosphate"                            ;   # 14
        10.0 ;       # "sedoheptulose-7-phosphate"                             ;   # 15

        10.0 ;       # "D-erythrose-4-phosphate"                               ;   # 16
        10.0 ;       # "D-fructose-1,6-bisphosphate"                           ;   # 17
        10.0 ;       # "dihydroxyacetone-phosphate"                            ;   # 18
        10.0 ;       # "3-phospho-D-glyceroyl-phosphate"                       ;   # 19
        10.0 ;       # "3-phospho-D-glycerate"                                 ;   # 20

        10.0 ;       # "D-glycerate-2-phosphate"                               ;   # 21
        10.0 ;       # "phosphoenolpyruvate"                                   ;   # 22
        10.0 ;       # "pyruvate"                                              ;   # 23
        10.0 ;       # "D-lactate"                                             ;   # 24
        25.0 ;       # "inorganic-phosphate"                                   ;   # 25
        10.0 ;       # "amp"                                                   ;   # 26
        100.0 ;       # "water"                                                 ;   # 27
        10. ;       # glyceraldehyde-3-phosphate                              ;   # 28

        Inf ;       # l-C
        Inf ;       # l-H
        Inf ;       # l-O
        Inf ;       # l-N
        Inf ;       # l-P
    ]

    # call main -
    (error_archive,element_error,parameter_archive) = main(species_list,
        species_dictionary,
        initial_composition_dictionary,
        lower_bound_array,
        upper_bound_array,
        initial_multiplier_dictionary,
        system_temperature_in_kelvin)

    # find the min error -
    min_index = indmin(error_archive)
    best = parameter_archive[:,min_index]
    writedlm("temp.out",best)

    push!(best_error_archive,error_archive[min_index])
    best_soln_archive = [best_soln_archive best]

    # grab the element error -
    push!(element_error_archive,element_error)

    # let the user know what is going on ...
    msg = "Completed $(run_index) of $(number_of_runs) trials ...\n"
    print(msg)

end

# dump to disk -
writedlm("best_error_archive.dat",best_error_archive)

# convert data back to mmol/L
for soln_index = 1:number_of_runs
    best_soln_archive[1:number_of_species,soln_index+1] = (FACTOR/VOLUME)*best_soln_archive[1:number_of_species,soln_index+1]
end

# dump the *corrected* (mmol/L) values to disk -
writedlm("best_soln_archive.dat",best_soln_archive[:,2:end])

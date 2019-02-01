include("Includes.jl")

# hard constants -
R = 8.314*(1/1000)                  # units: kJ mol^-1 K^-1
T = 298.15                          # units: K
penalty_weight = 10000000           # units: dimensionless

function calculate_penalty_term(mol_array, parameter_dictionary)

    # initialize -
    AM = parameter_dictionary["atom_matrix"]
    aV = parameter_dictionary["initial_atom_array"]

    # compute -
    error = AM*mol_array .- aV

    # compute penalty -
    number_of_constraints = length(error)
    penalty_array = zeros(number_of_constraints)
    for constraint_index = 1:number_of_constraints
        local_penalty_term = maximum([0,-1*error[constraint_index]])^2
        penalty_array[constraint_index] = local_penalty_term
    end

    # return -
    return penalty_array
end

function calculate_activity_array(mol_array, parameter_dictionary)

    # calculate the mol_total -
    mol_total = sum(mol_array)

    # compute the activity array -
    activity_array = Float64[]
    for mol_value in mol_array

        # compute -
        activity_value = mol_value/mol_total

        # cache -
        push!(activity_array,real(log(Complex(activity_value))))
    end

    # return -
    return activity_array
end

function calculate_gibbs_formation_term(mol_array, parameter_dictionary)

    # get the delta_f_array from the parameter_dictionary -
    gibbs_f_energy_array = parameter_dictionary["gibbs_energy_array_in_kj_mol"]

    # compute the gibbs formation array -
    gibbs_formation_array = mol_array.*gibbs_f_energy_array

    # return -
    return sum(gibbs_formation_array)
end

function scaled_gibbs_energy_function(mol_array, parameter_dictionary)

    # compute the formation term -
    delta_f_term = calculate_gibbs_formation_term(mol_array,parameter_dictionary)

    # calculate the activity array -
    activity_array = calculate_activity_array(mol_array, parameter_dictionary)

    # calculate the penalty term -
    penalty_array = calculate_penalty_term(mol_array, parameter_dictionary)

    # compute the gibbs energy -
    total_gibbs_energy = delta_f_term + (R*T)*sum(mol_array.*activity_array)

    # compute the total penalty -
    total_penalty = penalty_weight*sum(penalty_array)

    @show sum(penalty_array)

    # return -
    return total_gibbs_energy+total_penalty
end

function setup()

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
        "alpha-D-ribose-5-phosphate"                            ;   # 13
        "sedoheptulose-7-phosphate"                             ;   # 14
        "D-erythrose-4-phosphate"                               ;   # 15
        "D-fructose-1,6-bisphosphate"                           ;   # 16
        "dihydroxyacetone-phosphate"                            ;   # 17
        "3-phospho-D-glyceroyl-phosphate"                       ;   # 18
        "3-phospho-D-glycerate"                                 ;   # 19
        "D-glycerate-2-phosphate"                               ;   # 20
        "phosphoenolpyruvate"                                   ;   # 21
        "pyruvate"                                              ;   # 22
        "D-lactate"                                             ;   # 23
        "inorganic-phosphate"                                   ;   # 24
        "amp"                                                   ;   # 25
        "water"                                                 ;   # 26
        "glyceraldehyde-3-phosphate"                            ;   # 27
        "acetyl-coa"                                            ;   # 28
        "nicotinamide-adenine-dinucleotide"                     ;   # 29
        "nicotinamide-adenine-dinucleotide-reduced"             ;   # 30
        "oxoglutarate"                                          ;   # 31
        "succinyl-CoA"                                          ;   # 32
        "succinate"                                             ;   # 33
        "fumarate"                                              ;   # 34
        "malate"                                                ;   # 35
        "oxaloacetate"                                          ;   # 36
    ]

    # setup the composition (mmol/L) -
    scale_factor = 1.0
    initial_composition_dictionary = Dict()
    initial_composition_dictionary["glucose"] = 1.8
    initial_composition_dictionary["atp"] = 1.179
    initial_composition_dictionary["adp"] = 0.367
    initial_composition_dictionary["glucose-6-phosphate"] =  0.041
    initial_composition_dictionary["fructose-6-phosphate"] = 0.016
    initial_composition_dictionary["6-phospho-D-gluconate"] =  0.026
    initial_composition_dictionary["D-ribulose-5-phosphate"] = 0.087
    initial_composition_dictionary["5-phospho-alpha-D-ribose-1-diphosphate"] = 0.026
    initial_composition_dictionary["D-xylulose-5-phosphate"] =  0.013
    initial_composition_dictionary["nicotinamide-adenine-dinucleotide-phosphate"] = 1.075
    initial_composition_dictionary["nicotinamide-adenine-dinucleotide-phosphate-reduced"] = 2.0
    initial_composition_dictionary["co2"] = 1.2
    initial_composition_dictionary["alpha-D-ribose-5-phosphate"] = 0.078
    initial_composition_dictionary["sedoheptulose-7-phosphate"] =  0.016
    initial_composition_dictionary["D-erythrose-4-phosphate"] = 0.064
    initial_composition_dictionary["D-fructose-1,6-bisphosphate"] = 0.146
    initial_composition_dictionary["dihydroxyacetone-phosphate"] = 0.01
    initial_composition_dictionary["3-phospho-D-glyceroyl-phosphate"] = 0.023
    initial_composition_dictionary["3-phospho-D-glycerate"] =  0.01
    initial_composition_dictionary["D-glycerate-2-phosphate"] = 0.084
    initial_composition_dictionary["phosphoenolpyruvate"] = 0.011
    initial_composition_dictionary["pyruvate"] = 6.005
    initial_composition_dictionary["inorganic-phosphate"] = 1.00
    initial_composition_dictionary["D-lactate"] =  7.330
    initial_composition_dictionary["amp"] = 0.0831
    initial_composition_dictionary["water"] =  0.01
    initial_composition_dictionary["glyceraldehyde-3-phosphate"] = 0.00405
    initial_composition_dictionary["acetyl-coa"] = 0.00405
    initial_composition_dictionary["nicotinamide-adenine-dinucleotide"]= 1.0
    initial_composition_dictionary["nicotinamide-adenine-dinucleotide-reduced"]= 1.0
    initial_composition_dictionary["oxoglutarate"]= 1.0
    initial_composition_dictionary["succinyl-CoA"]= 1.0
    initial_composition_dictionary["succinate"]= 2.432
    initial_composition_dictionary["fumarate"]= 1.0
    initial_composition_dictionary["malate"]= 3.508
    initial_composition_dictionary["oxaloacetate"]= 1.0

    # number of states?
    number_of_species = length(species_list)

    # build the energy of formation array -
    gibbs_energy_array_in_kj_mol = buildGibbsEnergyOfFormationArray("./data/Database.json",species_list)

    # Build the species dictionary -
    species_dictionary = buildSpeciesDictionary("./data/Database.json")
    atom_matrix = buildAtomMatrix("./data/Database.json",species_list)

    # calculate the initial atom array -
    tmp_mol_array::Array{Float64,1} = Float64[]
    lower_bound_array = Float64[]
    upper_bound_array = Float64[]
    for (index,species_key) in enumerate(species_list)

        value = initial_composition_dictionary[species_key]
        push!(tmp_mol_array,value)
    end

    # initial atom vector (A)
    A = atom_matrix*tmp_mol_array

    # initial composition_array -
    total_initial_mol = sum(tmp_mol_array)

    # Setup the bounds -
    lower_bound_array = [

        1.759   ;     # "glucose"                                               ;   # 1
        0.00001 ;     # "fructose-6-phosphate"                                  ;   # 2
        0.00001 ;     # "glucose-6-phosphate"                                   ;   # 3
        0.1179  ;     # "atp"                                                   ;   # 4
        0.0367  ;     # "adp"                                                   ;   # 5

        0.00001 ;     # "6-phospho-D-gluconate"                                 ;   # 6
        0.00001 ;     # "D-ribulose 5-phosphate"                                ;   # 7
        0.00001 ;     # "5-phospho-alpha-D-ribose-1-diphosphate"                ;   # 8
        0.00001 ;     # "D-xylulose-5-phosphate"                                ;   # 9
        1.0     ;     # "nicotinamide-adenine-dinucleotide-phosphate"           ;   # 10

        1.2     ;     # "nicotinamide-adenine-dinucleotide-phosphate-reduced"   ;   # 11
        0.0001  ;     # "co2"                                                   ;   # 12
        0.00001 ;     # "alpha-D-Ribose-5-phosphate"                            ;   # 13
        0.00001 ;     # "sedoheptulose-7-phosphate"                             ;   # 14

        0.00001 ;     # "D-erythrose-4-phosphate"                               ;   # 15
        0.00001 ;     # "D-fructose-1,6-bisphosphate"                           ;   # 16
        0.00001 ;     # "dihydroxyacetone-phosphate"                            ;   # 17
        0.00001 ;     # "3-phospho-D-glyceroyl-phosphate"                       ;   # 18
        0.00001 ;     # "3-phospho-D-glycerate"                                 ;   # 19

        0.00001 ;     # "D-glycerate-2-phosphate"                               ;   # 20
        0.00001 ;     # "phosphoenolpyruvate"                                   ;   # 21
        0.6005  ;     # "pyruvate"                                              ;   # 22
        0.733   ;     # "D-lactate"                                             ;   # 23
        0.0001  ;     # "inorganic-phosphate"                                   ;   # 24
        0.00831 ;     # "amp"                                                   ;   # 25
        0.00001 ;     # "water"                                                 ;   # 26
        0.00001 ;     # "glyceraldehyde-3-phosphate"                            ;   # 27
        0.00001 ;     # acetyl-coa                                              ;   # 28

        0.1     ;     # nicotinamide-adenine-dinucleotide
        0.1     ;     # nicotinamide-adenine-dinucleotide-reduced
        0.1     ;     # oxoglutarate
        0.1     ;     # succinyl-CoA
        0.2432  ;     # succinate
        0.1     ;     # fumarate
        0.3508  ;     # malate
        0.1     ;     # oxaloacetate
    ]

    upper_bound_array = [

        2.00   ;       # "glucose"                                               ;   # 1
        10.0    ;       # "fructose-6-phosphate"                                  ;   # 2
        10.0    ;       # "glucose-6-phosphate"                                   ;   # 3
        11.79   ;       # "atp"                                                   ;   # 4
        3.674   ;       # "adp"                                                   ;   # 5
        10.0    ;       # "6-phospho-D-gluconate"                                 ;   # 6
        10.0    ;       # "D-ribulose-5-phosphate"                                ;   # 7
        10.0    ;       # "5-phospho-alpha-D-ribose-1-diphosphate"                ;   # 8
        10.0    ;       # "D-xylulose-5-phosphate"                                ;   # 9
        150.0   ;       # "nicotinamide-adenine-dinucleotide-phosphate"           ;   # 10

        150.0   ;       # "nicotinamide-adenine-dinucleotide-phosphate-reduced"   ;   # 11
        25.0    ;       # "co2"                                                   ;   # 12
        10.0    ;       # "alpha-D-Ribose-5-phosphate"                            ;   # 13
        10.0    ;       # "sedoheptulose-7-phosphate"                             ;   # 14

        10.0    ;       # "D-erythrose-4-phosphate"                               ;   # 15
        10.0    ;       # "D-fructose-1,6-bisphosphate"                           ;   # 16
        10.0    ;       # "dihydroxyacetone-phosphate"                            ;   # 17
        10.0    ;       # "3-phospho-D-glyceroyl-phosphate"                       ;   # 18
        10.0    ;       # "3-phospho-D-glycerate"                                 ;   # 19

        10.0    ;       # "D-glycerate-2-phosphate"                               ;   # 20
        10.0    ;       # "phosphoenolpyruvate"                                   ;   # 21
        60.05   ;       # "pyruvate"                                              ;   # 22
        73.295  ;       # "D-lactate"                                             ;   # 23
        25.0    ;       # "inorganic-phosphate"                                   ;   # 24
        0.8315  ;       # "amp"                                                   ;   # 25
        100.0   ;       # "water"                                                 ;   # 26
        10.0    ;       # glyceraldehyde-3-phosphate                              ;   # 27
        1.0     ;       # acetyl-coa                                              ;   # 28

        100     ;       # nicotinamide-adenine-dinucleotide
        100     ;       # nicotinamide-adenine-dinucleotide-reduced
        100     ;       # oxoglutarate
        100     ;       # succinyl-CoA
        24.32   ;       # succinate
        100     ;       # fumarate
        35.07   ;       # malate
        100     ;       # oxaloacetate

    ]

    # ======== DO NOT EDIT BELOW THIS LINE ====================================================== %
    parameter_dictionary = Dict()
    parameter_dictionary["gibbs_energy_array_in_kj_mol"] = gibbs_energy_array_in_kj_mol
    parameter_dictionary["species_dictionary"] = species_dictionary
    parameter_dictionary["number_of_species"] = number_of_species
    parameter_dictionary["initial_mol_array"] = tmp_mol_array
    parameter_dictionary["total_initial_number_of_mols"] = total_initial_mol
    parameter_dictionary["atom_matrix"] = atom_matrix
    parameter_dictionary["initial_atom_array"] = A
    parameter_dictionary["lower_bound_array"] = lower_bound_array
    parameter_dictionary["upper_bound_array"] = upper_bound_array
    return parameter_dictionary
    # =========================================================================================== %
end

function main()

    # build the parameter dictionary -
    parameter_dictionary = setup()

    # get a few things from the parameter_dictionary -
    number_of_species = parameter_dictionary["number_of_species"]
    initial_mol_array = parameter_dictionary["initial_mol_array"]
    initial_mol_total = parameter_dictionary["total_initial_number_of_mols"]

    # setup the call to the optimizer --
    objective_function(x) = scaled_gibbs_energy_function(x,parameter_dictionary)

    # get the lower, and upper_bounds -
    lower_bound_array = parameter_dictionary["lower_bound_array"]
    upper_bound_array = parameter_dictionary["upper_bound_array"]

    # make a call to the optim package -
    result = optimize(objective_function,lower_bound_array,upper_bound_array,initial_mol_array,Fminbox(LBFGS()))

    # get the min energy composition array -
    min_energy_composition_array = Optim.minimizer(result)

    # return -
    return min_energy_composition_array
end

# solve -
composition_array = main()

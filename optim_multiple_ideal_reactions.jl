include("Includes.jl")

# setup some global constants -
R_constant = 8.314 # units: J mol^-1 K^-1
system_temperature_in_kelvin = 298.15 # units: J mol^-1 K^-1

# helper functions -
function calculate_activity_array(extent_array,parameter_dictionary)

    # get stuff from the parameter_dictionary -
    stoichiometric_matrix = parameter_dictionary["stoichiometric_matrix"]
    number_of_reactions = parameter_dictionary["number_of_reactions"]
    number_of_species = parameter_dictionary["number_of_species"]

    # compute the scaled composition array -
    scaled_composition_array = calculate_scaled_composistion_array(extent_array,parameter_dictionary)

    # First, compute total number of mol at this extent -
    total_stcoeff_array = zeros(number_of_reactions)
    for reaction_index = 1:number_of_reactions
        stm_row_sum = sum(stoichiometric_matrix[:,reaction_index])
        total_stcoeff_array[reaction_index] = stm_row_sum
    end

    # total mol -
    current_total_mol = 1+sum(total_stcoeff_array.*extent_array)

    # compute activity -
    activity_array = zeros(number_of_species)
    for species_index = 1:number_of_species

        # activity value -
        activity_value = (1/current_total_mol)*scaled_composition_array[species_index]
        activity_array[species_index] = real(log(Complex(activity_value)))
    end

    # return -
    return activity_array
end

function calculate_scaled_composistion_array(extent_array, parameter_dictionary)

    # get stuff from the parameter_dictionary -
    initial_composition_array = parameter_dictionary["initial_composition_array"]
    stoichiometric_matrix = parameter_dictionary["stoichiometric_matrix"]

    # what is my initial mol total?
    initial_mol_total = sum(initial_composition_array[:,2])
    initial_mol_array = initial_composition_array[:,2]

    # compute the raw mol totals for this extent -
    raw_mol_array = initial_mol_array + stoichiometric_matrix*(extent_array*initial_mol_total)

    # scale by the total number of initail mols -
    scaled_mol_array = (1/initial_mol_total)*raw_mol_array

    # return -
    return scaled_mol_array
end

function calculate_gibbs_reaction_term(extent_array, parameter_dictionary)

    # get the scaled array of delta g's -
    scaled_delta_g_rxn_in_j_mol_array = parameter_dictionary["scaled_delta_g_rxn_in_j_mol_array"]

    # delta_g_term -
    delta_g_term = sum(extent_array.*scaled_delta_g_rxn_in_j_mol_array)

    # return -
    return delta_g_term
end

function scaled_gibbs_energy_function(extent_array, parameter_dictionary)

    # compute reaction term -
    delta_g_term = calculate_gibbs_reaction_term(extent_array,parameter_dictionary)

    # compute scaled scaled_composition_array -
    scaled_composition_array = calculate_scaled_composistion_array(extent_array, parameter_dictionary)

    # calculate the activity array -
    activity_array = calculate_activity_array(extent_array, parameter_dictionary)

    # compute the gibs energy -
    total_scaled_gibbs_energy = delta_g_term + sum(scaled_composition_array.*activity_array)

    # compute penalty -
    number_of_species = parameter_dictionary["number_of_species"]
    penalty_array = zeros(number_of_species)
    for species_index = 1:number_of_species
        local_penalty_term = maximum([0,-1*scaled_composition_array[species_index]])^2
        penalty_array[species_index] = local_penalty_term
    end

    # compute the total penalty -
    total_penalty = 1000*sum(penalty_array)

    # return -
    return total_scaled_gibbs_energy+total_penalty
end

function setup()

    # Build the species dictionary -
    species_dictionary = buildSpeciesDictionary("./data/TPI.json")

    # Build the reaction dictionary -
    reaction_dictionary = buildReactionDictionary("./data/TPI.json")

    # Setup species symbol array, setup the ICs -
    species_ic_array = [

        "glucose"               10.0                    ;   # 1
        "atp"                   20.0                    ;   # 2
        "glucose-6-phosphate"   0.1                     ;   # 3
        "adp"                   0.1                     ;   # 4
        "fructose-6-phosphate"  0.1                     ;   # 5
        "D-fructose-1,6-bisphosphate"   0.1             ;   # 6
        "glyceraldehyde-3-phosphate"    0.0             ;   # 7
        "dihydroxyacetone-phosphate"    10.0            ;   # 8
        "nicotinamide-adenine-dinucleotide" 10.0        ;   # 9
        "nicotinamide-adenine-dinucleotide-reduced" 0.1 ;   # 10
        "inorganic-phosphate" 10.0                      ;   # 11
        "3-phospho-D-glyceroyl-phosphate"   0.1         ;   # 12
        "3-phospho-D-glycerate" 0.1                     ;   # 13
        "D-glycerate-2-phosphate"    0.1                ;   # 14
        "phosphoenolpyruvate"   0.1                     ;   # 15
        "pyruvate"  0.1                                 ;   # 16
    ]

    species_symbol_array = species_ic_array[:,1]
    initial_mol_total = sum(species_ic_array[:,2])*(1/1000) # convert from mmol to mol
    initial_mol_array = species_ic_array[:,2]*(1/1000)  # convert from mmol to mol
    number_of_species = length(initial_mol_array)

    # order lookup dictionary -
    order_lookup_dictionary = Dict{String,Int}()
    for species_index = 1:number_of_species

        species_symbol = species_symbol_array[species_index]
        order_lookup_dictionary[species_symbol] = species_index
    end


    # What reactions id are we going to look at?
    reaction_key_array = String[]
    for (reaction_id,value) in reaction_dictionary
        push!(reaction_key_array,reaction_id)
    end
    number_of_reactions = length(reaction_key_array)

    # reaction_id = "2f8611d4-e3c8-4a89-bea9-f613f5574195"    # glucose + atp = glucose-6-phosphate + adp
    # push!(reaction_key_array,reaction_id)
    # reaction_id = "86b628c4-2eb3-43be-9014-9ac55f503503"    # glucose-6-phosphate = fructose-6-phosphate
    # push!(reaction_key_array,reaction_id)
    # reaction_id = "5ffef88f-98b1-4918-8188-cb4ac0e9f31e"
    # push!(reaction_key_array,reaction_id)
    # reaction_id = "5ffef88f-98b1-4918-8188-cb4ac0e9f31d"
    # push!(reaction_key_array,reaction_id)
    # reaction_id = "5ffef88f-98b1-4918-8188-cb4ac0e9f31c"
    # push!(reaction_key_array,reaction_id)

    # calculate the delta_g_rxn array -
    scaled_delta_g_rxn_in_j_mol_array = Float64[]
    for reaction_id in reaction_key_array

        # get my reaction and the associated species -
        reaction_object = reaction_dictionary[reaction_id]
        species_array = reaction_object.species_array

        # Build the stoichiometric dictionary -
        delta_g_rxn_in_j_mol = 0
        for species_reference_object::CRESpeciesReference in species_array

            # what species am I looking at?
            species_symbol = species_reference_object.symbol
            stoichiometric_coefficient = species_reference_object.stoichiometric_coefficient

            # lookup the dG of formation -
            dG_formation_in_j_mol = 1000*species_dictionary[species_symbol].delta_gibbs_in_kj_mol

            # compute the delta_g_rxn -
            delta_g_rxn_in_j_mol = delta_g_rxn_in_j_mol + stoichiometric_coefficient*dG_formation_in_j_mol
        end

        # scale by RT -
        scaled_delta_gibbs_reaction = (1/(R_constant*system_temperature_in_kelvin))*delta_g_rxn_in_j_mol
        push!(scaled_delta_g_rxn_in_j_mol_array,scaled_delta_gibbs_reaction)
    end

    # calculate the maximum extent array (upper_bound)
    upper_bound_array = ones(length(reaction_key_array))
    lower_bound_array = zeros(length(reaction_key_array))

    # compute the st matrix -
    stoichiometric_matrix = zeros(number_of_species,number_of_reactions)
    reaction_index = 1
    for reaction_key in reaction_key_array

        # get reaction object -
        reaction_object = reaction_dictionary[reaction_key]
        species_array = reaction_object.species_array
        for species_reference_object::CRESpeciesReference in species_array

            # what species am I looking at?
            species_symbol = species_reference_object.symbol
            stoichiometric_coefficient = species_reference_object.stoichiometric_coefficient

            # what species index are we?
            species_index = order_lookup_dictionary[species_symbol]

            # update stm value -
            stoichiometric_matrix[species_index,reaction_index] = stoichiometric_coefficient
        end

        # update reaction index -
        reaction_index = reaction_index + 1
    end


    # ======== DO NOT EDIT BELOW THIS LINE ======================================================= %
    parameter_dictionary = Dict()
    parameter_dictionary["initial_composition_array"] = species_ic_array
    parameter_dictionary["reaction_key_array"] = reaction_key_array
    parameter_dictionary["species_dictionary"] = species_dictionary
    parameter_dictionary["reaction_dictionary"] = reaction_dictionary
    parameter_dictionary["scaled_delta_g_rxn_in_j_mol_array"] = scaled_delta_g_rxn_in_j_mol_array
    parameter_dictionary["upper_bound_array"] = upper_bound_array
    parameter_dictionary["lower_bound_array"] = lower_bound_array
    parameter_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
    parameter_dictionary["total_initial_number_of_mols"] = initial_mol_total
    parameter_dictionary["initial_mol_array"] = initial_mol_array
    parameter_dictionary["number_of_reactions"] = number_of_reactions
    parameter_dictionary["number_of_species"] = number_of_species
    parameter_dictionary["species_symbol_array"] = species_symbol_array
    return parameter_dictionary
    # =========================================================================================== %
end

function main()

    # build the parameter dictionary -
    parameter_dictionary = setup()

    # get a few things from the parameter_dictionary -
    number_of_reactions = parameter_dictionary["number_of_reactions"]
    initial_mol_array = parameter_dictionary["initial_mol_array"]
    initial_mol_total = parameter_dictionary["total_initial_number_of_mols"]

    # setup the call to the optimizer --
    objective_function(x) = scaled_gibbs_energy_function(x,parameter_dictionary)

    # get the lower, and upper_bounds -
    lower_bound_array = parameter_dictionary["lower_bound_array"]
    upper_bound_array = parameter_dictionary["upper_bound_array"]

    # make a call to the optim package -
    result = optimize(objective_function,lower_bound_array,upper_bound_array,0.1*ones(number_of_reactions),Fminbox(LBFGS()))
    #result = optimize(objective_function,0,1)

    # get the extent -
    scaled_eq_extent_of_reaction = Optim.minimizer(result)

    # compute the final composition -
    final_scaled_composistion_array = calculate_scaled_composistion_array(scaled_eq_extent_of_reaction,parameter_dictionary)

    # return -
    return (scaled_eq_extent_of_reaction*initial_mol_total*1000,(final_scaled_composistion_array*initial_mol_total*1000), parameter_dictionary)
    # return parameter_dictionary
end

(extent,nf,pd) = main();

# create a results table -
number_of_species = pd["number_of_species"]
species_symbol_array = pd["species_symbol_array"]
initial_mol_array = pd["initial_mol_array"]*(1000)

results_table = Any[]
for species_index = 1:number_of_species

    species_symbol = species_symbol_array[species_index]
    record = [species_symbol  initial_mol_array[species_index] nf[species_index]]
    @show record

    push!(results_table,record)
    # results_table[species_index,1] = species_symbol
    # results_table[species_index,2] = initial_mol_array[species_index]
    # results_table[species_index,3] = nf[species_index]
end

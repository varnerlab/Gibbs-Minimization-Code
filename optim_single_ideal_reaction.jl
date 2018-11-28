include("Includes.jl")

function test_function(x,p)
    pval = p["value"]
    (pval - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

function scaled_gibbs_energy_function(extent_value,parameter_dictionary)

    # get stuff from parameter_dictionary -
    reaction_object = parameter_dictionary["reaction_object"]
    species_array = reaction_object.species_array
    stoichiometric_array = parameter_dictionary["stoichiometric_array"]
    initial_number_of_mol_dictionary = parameter_dictionary["initial_number_of_mol_dictionary"]
    initial_mol_total = parameter_dictionary["initial_mol_total"]
    scaled_delta_gibbs_reaction = parameter_dictionary["scaled_delta_gibbs_reaction"]

    # compute the total number of mols -
    species_mol_array = Float64[]
    for species_reference_object::CRESpeciesReference in species_array

        # what species am I looking at?
        species_symbol = species_reference_object.symbol

        # what is my stoichiometric_coefficient?
        stoichiometric_coefficient = stoichiometric_array[species_symbol]

        # what are the moles?
        mol_species_i = (initial_number_of_mol_dictionary[species_symbol]/initial_mol_total)+stoichiometric_coefficient*(extent_value[1])

        # cache this mol value -
        push!(species_mol_array,mol_species_i)
    end

    # ok, so what is the this total number of mols for this extent?
    total_number_of_mols = sum(species_mol_array)

    # compute the summation term -
    summation_term = 0
    species_index = 1
    for species_reference_object::CRESpeciesReference in species_array

        # get the dG of formation -
        species_symbol = species_reference_object.symbol

        # get the mol count for this species -
        mol_species_i = species_mol_array[species_index]

        # compute the summation term -
        summation_term = summation_term+mol_species_i*(log((mol_species_i/total_number_of_mols)))

        # update the counter -
        species_index = species_index + 1
    end

    # compute the gibbs total -
    gibbs_total = extent_value[1]*(scaled_delta_gibbs_reaction)+summation_term

    # return -
    return gibbs_total
end

function main(reaction_id,system_temperature_in_kelvin,species_dictionary,reaction_dictionary,initial_number_of_mol_dictionary)

    # constants -
    R_constant = 8.314

    # initialize -
    gibbs_array::Array{Float64,1} = Float64[]
    stoichiometric_array = Dict()

    # get my reaction and the associated species -
    reaction_object = reaction_dictionary[reaction_id]
    species_array = reaction_object.species_array

    # Build the stoichiometric dictionary -
    delta_g_rxn_in_j_mol = 0
    for species_reference_object::CRESpeciesReference in species_array

        # what species am I looking at?
        species_symbol = species_reference_object.symbol
        stoichiometric_coefficient = species_reference_object.stoichiometric_coefficient

        # add -
        stoichiometric_array[species_symbol] = stoichiometric_coefficient

        # lookup the dG of formation -
        dG_formation_in_j_mol = 1000*species_dictionary[species_symbol].delta_gibbs_in_kj_mol

        # compute the delta_g_rxn -
        delta_g_rxn_in_j_mol = delta_g_rxn_in_j_mol + stoichiometric_coefficient*dG_formation_in_j_mol
    end

    # scale by RT -
    scaled_delta_gibbs_reaction = (1/(R_constant*system_temperature_in_kelvin))*delta_g_rxn_in_j_mol

    # what is my initial mols total?
    initial_mol_total = 0.0
    for species_reference_object::CRESpeciesReference in species_array

        # what species am I looking at?
        species_symbol = species_reference_object.symbol

        # what is my initial ammount?
        initial_mol_number = initial_number_of_mol_dictionary[species_symbol]
        initial_mol_total = initial_mol_total + initial_mol_number
    end

    # find the max scaled extent -
    max_possible_extent_array = Float64[]
    for species_reference_object::CRESpeciesReference in species_array

        # what species am I looking at?
        species_symbol = species_reference_object.symbol

        if (species_reference_object.role == :reactant)

            # what is my stoichiometric_coefficient?
            stoichiometric_coefficient = stoichiometric_array[species_symbol]

            # calc max extent -
            max_extent = -(1/stoichiometric_coefficient)*(initial_number_of_mol_dictionary[species_symbol]/initial_mol_total)

            # cache -
            push!(max_possible_extent_array,max_extent)
        end
    end

    # setup the parameter dictionary -
    parameter_dictionary = Dict()
    parameter_dictionary["scaled_delta_gibbs_reaction"] = scaled_delta_gibbs_reaction
    parameter_dictionary["initial_mol_total"] = initial_mol_total
    parameter_dictionary["reaction_object"] = reaction_object
    parameter_dictionary["stoichiometric_array"] = stoichiometric_array
    parameter_dictionary["initial_number_of_mol_dictionary"] = initial_number_of_mol_dictionary


    # what is the lower and upper bounds on the extent?
    lower_bound = 0
    upper_bound = minimum(max_possible_extent_array)

    # setup the objective function -
    objective_function(x) = scaled_gibbs_energy_function(x,parameter_dictionary)

    # make a call to the optim package -
    result = optimize(objective_function,lower_bound,upper_bound)

    # get the results -
    scaled_eq_extent_of_reaction = Optim.minimizer(result)
    scaled_gibbs_energy = Optim.minimum(result)

    # convert back to normal units -
    eq_extent_of_reaction = initial_mol_total*scaled_eq_extent_of_reaction

    # gibbs -
    final_mol_array = Float64[]
    for species_reference_object::CRESpeciesReference in species_array

        # what species am I looking at?
        species_symbol = species_reference_object.symbol

        # what is my stoichiometric_coefficient?
        stoichiometric_coefficient = stoichiometric_array[species_symbol]

        # what is my initial ammount?
        initial_mol_number = initial_number_of_mol_dictionary[species_symbol]

        # mol -
        mol_value = initial_mol_number + (stoichiometric_coefficient*eq_extent_of_reaction)

        # cache -
        push!(final_mol_array,mol_value)
    end

    # final mol total -
    final_mol_total = sum(final_mol_array)
    mol_frac_array = (1/final_mol_total).*final_mol_array

    species_index = 1
    Q = 1
    for species_reference_object::CRESpeciesReference in species_array

        # what species am I looking at?
        species_symbol = species_reference_object.symbol

        # what is my stoichiometric_coefficient?
        stoichiometric_coefficient = stoichiometric_array[species_symbol]
        Q = Q*(mol_frac_array[species_index])^(stoichiometric_coefficient)
        species_index = species_index + 1
    end

    @show delta_g_rxn_in_j_mol,(R_constant*system_temperature_in_kelvin)*log(Q)

    abs_gibbs_energy = delta_g_rxn_in_j_mol + (R_constant*system_temperature_in_kelvin)*log(Q)

    # return -
    return (eq_extent_of_reaction,abs_gibbs_energy)
end


# Build the species dictionary -
species_dictionary = buildSpeciesDictionary("./data/Database.json")

# Build the reaction dictionary -
reaction_dictionary = buildReactionDictionary("./data/Database.json")

# What reaction id are we going to look at?
reaction_id = "86b628c4-2eb3-43be-9014-9ac55f503503"
initial_number_of_mol_dictionary = Dict()
initial_number_of_mol_dictionary["glucose-6-phosphate"] = 10.0
initial_number_of_mol_dictionary["fructose-6-phosphate"] = 1.0

# reaction_id = "5ffef88f-98b1-4918-8188-cb4ac0e9f31b"
# initial_number_of_mol_dictionary = Dict()
# initial_number_of_mol_dictionary["citrate"] = 1.0
# initial_number_of_mol_dictionary["isocitrate"] = 1e-9

# reaction_id = "5ffef88f-98b1-4918-8188-cb4ac0e9f31c"
# initial_number_of_mol_dictionary = Dict()
# initial_number_of_mol_dictionary["dihydroxyacetone-phosphate"] = 10.0
# initial_number_of_mol_dictionary["glyceraldehyde-3-phosphate"] = 1.0

# reaction_id = "5ffef88f-98b1-4918-8188-cb4ac0e9f31d"
# initial_number_of_mol_dictionary = Dict()
# initial_number_of_mol_dictionary["D-fructose-1,6-bisphosphate"] = 10.0
# initial_number_of_mol_dictionary["dihydroxyacetone-phosphate"] = 1.0
# initial_number_of_mol_dictionary["glyceraldehyde-3-phosphate"] = 1.0

# What is the system T (in K)
system_temperature_in_kelvin = 298.15

# calculate the gibbs as a function of e -
(extent_of_reaction,gibbs_energy) = main(reaction_id,system_temperature_in_kelvin,species_dictionary,reaction_dictionary,initial_number_of_mol_dictionary)

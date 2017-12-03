# include -
include("Includes.jl")


function main(reaction_id,system_temperature_in_kelvin,species_dictionary,reaction_dictionary,initial_number_of_mol_dictionary)

    # constants -
    R_constant = 8.314

    # initialize -
    extent_of_reaction = 0.0
    K = 0.0
    stoichiometric_array = Dict()

    # get my reaction -
    reaction_object = reaction_dictionary[reaction_id]

    # calculate the dG for this reaction -
    delta_gibbs_reaction_in_j_mol = 0.0
    species_array = reaction_object.species_array
    for species_reference_object::CRESpeciesReference in species_array

        # what species am I looking at?
        species_symbol = species_reference_object.symbol
        stoichiometric_coefficient = species_reference_object.stoichiometric_coefficient

        # add -
        stoichiometric_array[species_symbol] = stoichiometric_coefficient

        # lookup the dG -
        dG_formation = 1000*species_dictionary[species_symbol].delta_gibbs_in_kj_mol

        # calculate dG_reaction -
        delta_gibbs_reaction_in_j_mol = delta_gibbs_reaction_in_j_mol+stoichiometric_coefficient*dG_formation
    end

    # scale by RT -
    scaled_delta_gibbs_reaction = (1/(R_constant*system_temperature_in_kelvin))*delta_gibbs_reaction_in_j_mol

    # What is my K?
    K = exp(-1*scaled_delta_gibbs_reaction)

    # now estimate the eq extent of reaction -
    extent_of_reaction_array = collect(linspace(0,0.999,1000))
    error_array::Array{Float64,1} = Float64[]
    for extent_value in extent_of_reaction_array

        # initialize -
        final_state_array = Dict()

        # what is my initial mols total?
        initial_mol_total = 0.0
        for species_reference_object::CRESpeciesReference in species_array

            # what species am I looking at?
            species_symbol = species_reference_object.symbol

            # what is my initial ammount?
            initial_mol_number = initial_number_of_mol_dictionary[species_symbol]
            initial_mol_total = initial_mol_total + initial_mol_number
        end

        # calculate the state -
        mol_total = 0
        for species_reference_object::CRESpeciesReference in species_array

            # what species am I looking at?
            species_symbol = species_reference_object.symbol

            # what is my initial ammount?
            initial_mol_number = initial_number_of_mol_dictionary[species_symbol]

            # what is my stoichiometric_coefficient?
            stoichiometric_coefficient = stoichiometric_array[species_symbol]

            # what is my final state?
            value = initial_mol_number+stoichiometric_coefficient*(extent_value*initial_mol_total)
            final_state_array[species_symbol] = value

            # update total -
            mol_total = mol_total+value
        end

        # calculate the equilibrium constant -
        K_computed = 1
        for species_reference_object::CRESpeciesReference in species_array

            # what species am I looking at?
            species_symbol = species_reference_object.symbol
            stoichiometric_coefficient = stoichiometric_array[species_symbol]

            # compute the mol fraction -
            mol_fraction_species_i = initial_number_of_mol_dictionary[species_symbol]/(initial_mol_total)+stoichiometric_coefficient*extent_value

            # update the K_computed -
            K_computed = K_computed*(mol_fraction_species_i)^stoichiometric_coefficient
        end

        # what is the diff between K and K_computed?
        error = K - K_computed
        push!(error_array,error)
    end

    # sort -
    idx_zero = find(error_array.>0)
    eq_extent = extent_of_reaction_array[idx_zero[end]]

    return (error_array,extent_of_reaction_array,K,eq_extent)
end

# execute the main -

# Build the species dictionary -
species_dictionary = buildSpeciesDictionary("./data/Database.json")

# Build the reaction dictionary -
reaction_dictionary = buildReactionDictionary("./data/Database.json")

# What reaction id are we going to look at?
reaction_id = "86b628c4-2eb3-43be-9014-9ac55f503503"
initial_number_of_mol_dictionary = Dict()
initial_number_of_mol_dictionary["glucose-6-phosphate"] = 10.0
initial_number_of_mol_dictionary["fructose-6-phosphate"] = 1e-9

# What reaction id are we going to look at?
# reaction_id = "2f8611d4-e3c8-4a89-bea9-f613f5574195"
# initial_number_of_mol_dictionary = Dict()
# initial_number_of_mol_dictionary["glucose"] = 10.0
# initial_number_of_mol_dictionary["atp"] = 10.0
# initial_number_of_mol_dictionary["glucose-6-phosphate"] = 1e-9
# initial_number_of_mol_dictionary["adp"] = 1e-9

# What is the system T (in K)
system_temperature_in_kelvin = 298.15

# calculate the equilibrium extent of conversion -
(error_array,extent_of_reaction_array,K,eq_extent) = main(reaction_id,system_temperature_in_kelvin,species_dictionary,reaction_dictionary,initial_number_of_mol_dictionary)

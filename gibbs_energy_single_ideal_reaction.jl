# include -
include("Includes.jl")

function main(reaction_id,system_temperature_in_kelvin,species_dictionary,reaction_dictionary,initial_number_of_mol_dictionary)

    # constants -
    R_constant = 8.314

    # initialize -
    gibbs_array::Array{Float64,1} = Float64[]
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



    upper_bound = minimum(max_possible_extent_array)
    extent_of_reaction_array = collect(linspace(0,upper_bound,1000))

    # main loop -
    for extent_value in extent_of_reaction_array

        # initialize -
        final_state_array = Dict()

        summation_term = 0
        for species_reference_object::CRESpeciesReference in species_array

            # what species am I looking at?
            species_symbol = species_reference_object.symbol

            # what is my stoichiometric_coefficient?
            stoichiometric_coefficient = stoichiometric_array[species_symbol]

            # what are the moles?
            mol_species_i = (initial_number_of_mol_dictionary[species_symbol]/initial_mol_total)+stoichiometric_coefficient*(extent_value)

            @show (species_symbol,mol_species_i,extent_value)

            # compute the summation term -
            summation_term = summation_term+mol_species_i*log(e,mol_species_i)
        end

        # finally ... calculate the total gibbs energy -
        gibbs_total = scaled_delta_gibbs_reaction*(extent_value)+summation_term
        push!(gibbs_array,gibbs_total)
    end

    return (extent_of_reaction_array,gibbs_array)
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
# initial_number_of_mol_dictionary["glucose"] = 1.0
# initial_number_of_mol_dictionary["atp"] = 1.0
# initial_number_of_mol_dictionary["glucose-6-phosphate"] = 1e-9
# initial_number_of_mol_dictionary["adp"] = 1e-9


# What is the system T (in K)
system_temperature_in_kelvin = 298.15

# calculate the gibbs as a function of e -
(extent_of_reaction_array,gibbs_array) = main(reaction_id,system_temperature_in_kelvin,species_dictionary,reaction_dictionary,initial_number_of_mol_dictionary)

# make a plot -
plot(extent_of_reaction_array,gibbs_array,"k",linewidth=2.0)

# axis labels -
xlabel_text = L"Scaled extent of reaction $\epsilon^{\prime}$ [AU]"
xlabel(xlabel_text,fontsize=16)
ylabel("Scaled Gibbs Energy [AU]",fontsize=16)

# find the min point -
idx_min = indmin(gibbs_array)
plot(extent_of_reaction_array[idx_min],gibbs_array[idx_min],"o",linewidth=4.0,mec="k",mfc="w",ms="10")

x_position = 0.01*maximum(extent_of_reaction_array)
y_position = 0.95*maximum(gibbs_array)

text_label = L"Min scaled extent of reaction $\epsilon^{\prime}$: "*string(extent_of_reaction_array[idx_min])
text(x_position,y_position,text_label,fontsize=16)

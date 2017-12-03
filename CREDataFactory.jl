function buildSpeciesDictionary(path_to_species_file::AbstractString)

    # Initialize species dictionary -
    species_dictionary::Dict{AbstractString,CRESpecies} = Dict()

    # load the database file w/the species information in it -
    raw_database_dict = JSON.parsefile(path_to_species_file)

    # ok, so we need to go though and build up the array for VLESpecies objects -
    raw_species_array = raw_database_dict["species_array"]
    for raw_species_data in raw_species_array

        # Build new species type -
        species_object = CRESpecies()

        # populate -
        species_object.symbol = raw_species_data["symbol"]
        species_object.delta_gibbs_in_kj_mol = parse(Float64,raw_species_data["delta_gibbs_in_kj_mol"])

        # get the element_array -
        species_element_dictionary::Dict{String,Float64} = Dict()
        element_array = raw_species_data["element_array"]
        for key in keys(element_array)

            # lookup the atom count -
            species_element_dictionary[key] = parse(Float64,element_array[key])
        end

        # grab the element dictionary for this species -
        species_object.element_dictionary = species_element_dictionary

        # store -
        species_dictionary[species_object.symbol] = species_object
    end

    return species_dictionary
end

function buildGibbsEnergyOfFormationArray(path_to_species_file::AbstractString,species_list::Array{String,1})

    # hard constants -
    VOLUME = 1.186e-13
    FACTOR = 1e-12

    # build the species_dictionary -
    species_dictionary = buildSpeciesDictionary(path_to_species_file)

    # initialize -
    gibbs_energy_array_in_j_mol::Array{Float64,1} = Float64[]

    # iterate through the species, and grab the dG of formation -
    for species_symbol in species_list

        # get species object -
        species_object = species_dictionary[species_symbol]

        # get dG -
        value = (FACTOR)*species_object.delta_gibbs_in_kj_mol
        push!(gibbs_energy_array_in_j_mol,value)

    end

    return gibbs_energy_array_in_j_mol
end

function buildAtomMatrix(path_to_species_file::AbstractString,species_list::Array{String,1})

    # load the database file w/the species information in it -
    raw_database_dict = JSON.parsefile(path_to_species_file)

    # build the species_dictionary -
    species_dictionary = buildSpeciesDictionary(path_to_species_file)

    # how many atoms?
    number_of_atoms = parse(Int64,raw_database_dict["number_of_atoms"])

    # how many species?
    number_of_species = length(species_list)

    # initialize -
    atom_matrix = zeros(number_of_atoms,number_of_species)

    # what is the atom order?
    atom_symbol_array = raw_database_dict["atom_symbol_array"]
    for (index_outer,atom_key) in enumerate(atom_symbol_array)

        for (index_inner,species_key) in enumerate(species_list)

            # Get the species object -
            species_object = species_dictionary[species_key]

            # what is my element_dictionary?
            element_dictionary = species_object.element_dictionary

            # fill in -
            atom_matrix[index_outer,index_inner] = element_dictionary[atom_key]
        end
    end

    return atom_matrix
end

function buildReactionDictionary(path_to_species_file::AbstractString)

    # Initialize species dictionary -
    reaction_dictionary::Dict{AbstractString,CREReaction} = Dict()

    # load the database file w/the species information in it -
    raw_database_dict = JSON.parsefile(path_to_species_file)

    # ok, so we need to go though and build up the array for VLESpecies objects -
    raw_reaction_array = raw_database_dict["reaction_array"]
    for raw_reaction_data in raw_reaction_array

        # Build new reaction type -
        reaction_object = CREReaction()

        # populate -
        reaction_object.symbol = raw_reaction_data["symbol"]
        reaction_object.id = raw_reaction_data["id"]

        # build the reactant array -
        reactant_array::Array{CRESpeciesReference,1} = CRESpeciesReference[]
        raw_reactant_array = raw_reaction_data["reactants"]
        for reactant in raw_reactant_array

            species_reference = CRESpeciesReference()
            species_reference.symbol = reactant["symbol"]
            species_reference.stoichiometric_coefficient = -1.0*parse(Float64,reactant["stoichiometric_coefficient"])
            species_reference.role = :reactant

            # cache -
            push!(reactant_array,species_reference)
        end

        # build the product array -
        raw_product_array = raw_reaction_data["products"]
        for product in raw_product_array

            species_reference = CRESpeciesReference()
            species_reference.symbol = product["symbol"]
            species_reference.stoichiometric_coefficient = parse(Float64,product["stoichiometric_coefficient"])
            species_reference.role = :product

            # cache -
            push!(reactant_array,species_reference)
        end

        reaction_object.species_array = reactant_array

        # store -
        reaction_dictionary[reaction_object.id] = reaction_object
    end

    return reaction_dictionary
end

mutable struct CRESpecies

    # species symbol -
    symbol::AbstractString
    delta_gibbs_in_kj_mol::Float64
    element_dictionary::Dict{String,Float64}

    function CRESpecies()
        this = new()
    end
end

mutable struct CRESpeciesReference

    # species symbol -
    symbol::AbstractString
    stoichiometric_coefficient::Float64
    role::Symbol #{:reactant,:product}

    function CRESpeciesReference()
        this = new()
    end
end

mutable struct CREReaction

    # species symbol -
    symbol::AbstractString
    id::AbstractString

    # reactant and product array -
    species_array::Array{CRESpeciesReference,1}

    function CREReaction()
        this = new()
    end
end

mutable struct CREProblemDataModel

    # dictionaries -
    species_dictionary::Dict{AbstractString,CRESpecies}
    initial_composition_dictionary::Dict{String,Float64}
    initial_multiplier_dictionary::Dict{String,Float64}
    initial_atom_dictionary::Dict{String,Float64}
    species_list::Array{String,1}
    system_temperature_in_kelvin::Float64

    atom_matrix::Array{Float64,2}
    total_atom_array::Array{Float64,1}
    gibbs_energy_array_in_j_mol::Array{Float64,1}

    parameter_lower_bound_array::Array{Float64}
    parameter_upper_bound_array::Array{Float64}

    number_of_species::Int

    function CREProblemDataModel()
        this = new()
    end
end

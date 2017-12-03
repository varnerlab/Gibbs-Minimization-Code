include("Includes.jl")

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
    VOLUME = 1.186e-13
    FACTOR = 1e-12

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
    size(atom_matrix)
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



function main(species_list,species_dictionary,initial_composition_dictionary,initial_multiplier_dictionary,system_temperature_in_kelvin)

    # hard constants -
    VOLUME = 1.186e-13
    FACTOR = 1e-12

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

    # Setup the bounds -
    LOWER_BOUND = [

        0.0 ;     # "glucose-6-phosphate"   ;   # 1
        0.0 ;     # "fructose-6-phosphate"  ;   # 2


        -Inf ;    # l-C
        -Inf ;    # l-H
        -Inf ;    # l-O
        -Inf ;    # l-N
        -Inf ;    # l-P
    ]

    UPPER_BOUND = [

        Inf ;    # "glucose-6-phosphate"   ;   # 1
        Inf ;    # "fructose-6-phosphate"  ;   # 2


        Inf ;     # l-C
        Inf ;     # l-H
        Inf ;     # l-O
        Inf ;     # l-N
        Inf ;     # l-P
    ]

    problem_data_model.parameter_lower_bound_array = (LOWER_BOUND)*(VOLUME/FACTOR)
    problem_data_model.parameter_upper_bound_array = (UPPER_BOUND)*(VOLUME/FACTOR)

    # call the minimizer -
    (error_archive,parameter_archive) = sa_estimate(objective_function,initial_parameter_guess,problem_data_model)

    # return -
    return (error_archive,parameter_archive)
end

# hard constants -
VOLUME = 1.186e-13
FACTOR = 1e-12

# specify the metabolites I want to optimize -
species_list = [
    "glucose-6-phosphate"                                   ;   # 1
    "fructose-6-phosphate"                                  ;   # 2
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
best_soln_archive = zeros(number_of_species+5,1)
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
        initial_composition_dictionary["glucose-6-phosphate"] =  1.00*(VOLUME/FACTOR)
        initial_composition_dictionary["fructose-6-phosphate"] = 0.0001*(VOLUME/FACTOR)

        # setup multiplier -
        initial_multiplier_dictionary["C"] = 100.0
        initial_multiplier_dictionary["H"] = 100.0
        initial_multiplier_dictionary["O"] = 100.0
        initial_multiplier_dictionary["N"] = 100.0
        initial_multiplier_dictionary["P"] = 100.0
    end

    # call main -
    (error_archive,parameter_archive) = main(species_list,species_dictionary,initial_composition_dictionary,initial_multiplier_dictionary,system_temperature_in_kelvin)

    # find the min error -
    min_index = indmin(error_archive)
    best = parameter_archive[:,min_index]
    writedlm("temp.out",best)

    push!(best_error_archive,error_archive[min_index])
    best_soln_archive = [best_soln_archive best]
end

# dump to disk -
writedlm("best_error_archive.dat",best_error_archive)
writedlm("best_soln_archive.dat",best_soln_archive[:,2:end])

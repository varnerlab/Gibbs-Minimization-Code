function sa_estimate(objective_function::Function,
    initial_parameter_guess::Array{Float64,1},
    problem_data_model::CREProblemDataModel)

    # How many parameters do we have?
    number_of_parameters = length(initial_parameter_guess)

    # setup some problem specific parameters -
    max_number_of_function_evaluations_per_temp = 1000
    max_number_of_temperature_adjustments = 500
    max_number_of_starts = 1
    should_we_continue = true
    iteration_counter = 1

    # storage -
    error_archive = zeros(1)
    parameter_archive = zeros(number_of_parameters)

    # calculate the initial error -
    best_error = objective_function(initial_parameter_guess,problem_data_model)
    best_parameter_array = initial_parameter_guess
    temperature = abs(log(e,best_error))

    # main loop -
    while (should_we_continue)

        for temperature_counter  = 1:max_number_of_temperature_adjustments
            for function_eval_counter = 1:max_number_of_function_evaluations_per_temp

                # generate a new parameter guess -
                new_parameter_guess = neighbor_function(best_parameter_array,problem_data_model,0.005)

                # Evaluate this parameter set -
                new_error = objective_function(new_parameter_guess,problem_data_model)

                # compute the error diff
                error_difference = (new_error - best_error)

                # compute the acceptance probality -
                acceptance_probability = acceptance_probability_function(error_difference,temperature)

                # should we keep this parameter set or not?
                if (acceptance_probability>rand())

                    # cache -
                    error_archive = [error_archive new_error]
                    parameter_archive = [parameter_archive new_parameter_guess]

                    # reset -
                    best_error = new_error
                    best_parameter_array = new_parameter_guess

                    msg = "Found new best error: $(best_error)\n"
                    print(msg)
                end
            end

            # update my temperature -
            temperature = cooling_function(temperature)
        end

        # we've gone through a run, should we go around again?
        if (iteration_counter >= max_number_of_starts)
            should_we_continue = false
        else

            # reset -
            best_parameter_array = neighbor_function(best_parameter_array,problem_data_model,0.5)
            best_error = objective_function(best_parameter_array,problem_data_model)
            temperature = abs(log(e,best_error))

            # update the counter -
            iteration_counter = iteration_counter + 1
        end
    end

    # TODO: We should check to see if we have empty arrays to correct

    return (error_archive[:,2:end],parameter_archive[:,2:end])
end

function acceptance_probability_function(error_difference,temperature)

    if (error_difference<0)
        return 1.0
    else
        return (exp(-error_difference/temperature))
    end
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.90
  return alpha*temperature
end


function neighbor_function(parameter_array,problem_data_model::CREProblemDataModel,SIGMA)

  number_of_parameters = length(parameter_array)

  # calculate new parameters -
  new_parameter_array = parameter_array.*(1+SIGMA*randn(number_of_parameters))

  # get the lower and upper bounds from the problem model -
  LOWER_BOUND = problem_data_model.parameter_lower_bound_array
  UPPER_BOUND = problem_data_model.parameter_upper_bound_array

  # return the corrected parameter arrays -
  return parameter_bounds_function(new_parameter_array,LOWER_BOUND,UPPER_BOUND)
end

function parameter_bounds_function(parameter_array,lower_bound_array,upper_bound_array)

  # reflection_factor -
  epsilon = 0.01

  # iterate through and fix the parameters -
  new_parameter_array = copy(parameter_array)
  for (index,value) in enumerate(parameter_array)

    lower_bound = lower_bound_array[index]
    upper_bound = upper_bound_array[index]

    if (value<lower_bound)
      new_parameter_array[index] = lower_bound
    elseif (value>upper_bound)
      new_parameter_array[index] = upper_bound
    end
  end

  return new_parameter_array
end

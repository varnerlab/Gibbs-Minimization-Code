include("Includes.jl")

function calculate_basis_null_space(matrix)

    # calculate the basis of matrix -
    # first, rref the matrix -
    reduced_matrix = rref(matrix)

    # next, which cols are constrained, and which are free?
    (number_of_rows,number_of_cols) = size(reduced_matrix)

    # how many free variables do we have?
    number_of_free_variables = 0
    free_variable_index_array = Int64[]
    bound_variable_index_array = Int64[]
    for col_index = 1:number_of_cols

        # look for cols w/a single non-zero element -
        idx_one = findall((reduced_matrix[:,col_index].== 1.0))
        idx_z = findall((reduced_matrix[:,col_index].== 0.0))
        if (length(idx_one) == 1 && length(idx_z) == (number_of_rows - 1))

            # capture the *free* index -
            push!(bound_variable_index_array,col_index)
        else

            # update the number of free variables -
            number_of_free_variables = number_of_free_variables + 1

            # capture the *bound* index -
            push!(free_variable_index_array,col_index)
        end
    end

    # ok, intialize the basis array -
    basis_array = zeros(number_of_cols,number_of_free_variables)

    # add the free variable indexes to the array -
    current_column_index = 1
    for free_variable_index in free_variable_index_array

        # fill in the free variables -
        basis_array[free_variable_index,current_column_index] = 1.0

        # next column -
        current_column_index = current_column_index + 1
    end


    # fill in other values -
    number_of_bound_variables = length(bound_variable_index_array)
    for row_index = 1:number_of_bound_variables

        # what is the current bound_variable_index?
        bound_variable_index = bound_variable_index_array[row_index]

        # fill in the missing values -
        for col_index = 1:number_of_free_variables
            free_variable_index = free_variable_index_array[col_index]
            basis_array[bound_variable_index,col_index] = -1*reduced_matrix[row_index,free_variable_index]
        end
    end

    return (bound_variable_index_array,free_variable_index_array,basis_array)
end

STM = [
    0 -0.5 -1 1 ;
    -1 -0.5 1 0 ;
    -1 -1 0 1 ;
]

# STM = [
#
#     -1 -1 1 1 0 0 ;
#     -1 0 0 1 -1 1 ;
#     0 -1 1 0 1 -1 ;
# ]

TSTM = Matrix(transpose(STM))

(bound_variable_index_array,free_variable_index_array,basis_array) = calculate_basis_null_space(TSTM)

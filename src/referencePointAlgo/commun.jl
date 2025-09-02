# =============================================================================
"""
    lex_optimize_knapsack(  solver::DataType,             
                            p::Matrix{Int64},   w::Vector{Int64},   c::Int64, 
                            obj_index::Int64
                            )

    Compute the lexicographic optimal solution for the bi-objective 01UKP
"""
function lex_optimize_knapsack( 
            solver::DataType,             
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64, 
            obj_index::Int64
        )

    d,n = size(p)
    model = Model(solver)
    @variable(model, x[1:n], Bin)

    # step 1 : minimize the primary objective ---------------------------------

    if obj_index == 1
        @objective(model, Max, sum(p[1,j] * x[j] for j in 1:n))
    elseif obj_index == 2
        @objective(model, Max, sum(p[2,j] * x[j] for j in 1:n))
    end
    @constraint(model, sum(w[j] * x[j] for j in 1:n) ≤ c)

    optimize!(model)
    @assert is_solved_and_feasible(model) "Error: optimal solution not found"

    f1_val = sum(p[1,j] * value(x[j]) for j in 1:n)
    f2_val = sum(p[2,j] * value(x[j]) for j in 1:n)

    # step 2 : minimize the secondary objective -------------------------------

    if obj_index == 1
        @objective(model, Max, sum(p[2,j] * x[j] for j in 1:n))
        @constraint(model, sum(p[1,j] * x[j] for j in 1:n) ≥ f1_val)
    else
        @objective(model, Max, sum(p[1,j] * x[j] for j in 1:n))
        @constraint(model, sum(p[2,j] * x[j] for j in 1:n) ≥ f2_val)
    end

    optimize!(model)
    @assert is_solved_and_feasible(model) "Error: optimal solution not found"

    f1_val = sum(p[1,j] * value(x[j]) for j in 1:n)
    f2_val = sum(p[2,j] * value(x[j]) for j in 1:n)
    
    return ( Int.(round.(f1_val)), Int.(round.(f2_val)) )
end


# =============================================================================
"""
    SolveAugmentedWeightedTchebycheff(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64,  
            λ::Vector{Float64}, rp::Vector{Int64}, ρ::Float64)

    Compute one nondominated point by application of the Tchebycheff theory.
    Application to the bi-objective 01 unidimensionnal knapsack problem. 
        
    In  Steuer, R.E., Choo, EU. An interactive weighted Tchebycheff procedure 
        for multiple objective programming. Mathematical Programming 26, 326–344 (1983).    
"""
function  SolveAugmentedWeightedTchebycheff(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64,  
            λ::Vector{Float64}, rp::Vector{Int64}, ρ::Float64
        )
    
    # number of objectives and number of variables        
    d,n = size(p)

    # define a JuMP model
    model = Model(solver)  

    # ------------------------------------------------------------------------------
    # Part related to the MO 01UKP model
    @variable(model, x[1:n], Bin)                                    # the n binary variables of 01UKP
    @constraint(model, sum(w[i] * x[i] for i in 1:n) ≤ c)            # the constraint of 01UKP

    # declare the objective functions
    @expression(model, z[k=1:d], sum(p[k,j] * x[j] for j in 1:n))    # the d objectives of 01UKP

    # ------------------------------------------------------------------------------
    # Part related to the Tchebycheff model
    @variable(model, α ≥ 0)        
    @objective(model, Min, α + ρ * sum((rp[k] -  z[k]) for k=1:d))

    # add the d constraints for a given λ to the model
    @constraint(model, con[k=1:d], α ≥  λ[k] * (rp[k] -  z[k]))

    set_silent(model)
    solve_time_sec = 0.0
    optimize!(model)
    @assert is_solved_and_feasible(model) "Error: optimal solution not found"
    solve_time_sec += solve_time(model)
    
    f1_val = sum(p[1,j] * value(x[j]) for j in 1:n)
    f2_val = sum(p[2,j] * value(x[j]) for j in 1:n)

    return round(Int, value(f1_val)), round(Int, value(f2_val)) #,x_sol = value.(x)
end

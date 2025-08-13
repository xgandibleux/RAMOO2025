using Printf
using JuMP, GLPK                         # for solving MILP (I)
#using HiGHS, Gurobi, CPLEX              # for solving MILP (II)
using Printf

#=
"""
Structure pour représenter un point dans l'espace objectif
"""
struct Point2D
    z1::Float64
    z2::Float64
end


Base.:(==)(p1::Point2D, p2::Point2D) = (p1.z1 == p2.z1) && (p1.z2 == p2.z2)
Base.isless(p1::Point2D, p2::Point2D) = p1.z1 < p2.z1 || (p1.z1 == p2.z1 && p1.z2 < p2.z2)

"""
Structure pour représenter un problème d'optimisation bicritère discret
"""
abstract type BicriteriaDiscreteProlem end

"""
Problème de sac à dos bicritère comme exemple d'implémentation
"""
mutable struct BicriteriaKnapsackProblem <: BicriteriaDiscreteProlem
    c1::Vector{Int}  # coefficients de coût pour objectif 1
    c2::Vector{Int}  # coefficients de coût pour objectif 2
    a::Vector{Int}   # coefficients de contrainte (profits)
    b::Int           # seuil minimum de profit
    n::Int           # nombre d'items
    
    function BicriteriaKnapsackProblem(c1::Vector{Int}, c2::Vector{Int}, a::Vector{Int}, b::Int)
        n = length(c1)
        @assert length(c2) == n && length(a) == n "Toutes les dimensions doivent être égales"
        new(c1, c2, a, b, n)
    end
end
=#

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


"""
    SolveAugmentedWeightedTchebycheff(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64,  
            λ::Vector{Float64}, rp::Vector{Int64}, ρ::Float64)

    Compute one nondominated point by application of the Tchebycheff theory
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


"""
    SteuerChoo(            
        solver::DataType,             
        p::Matrix{Int64},  w::Vector{Int64},  c::Int64
    )

    Interactive method proposed by Steuer and Choo in 1983
    (interactive part has to be implemented)

    In  Steuer, R.E., Choo, EU. An interactive weighted Tchebycheff procedure 
        for multiple objective programming. Mathematical Programming 26, 326–344 (1983).
"""
function SteuerChoo(            
            solver::DataType,             
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64
        )

    # Compute lexicographic optimal solutions z1 and z2 wrt. f1 and f2, respectively
    zOptf1 = lex_optimize_knapsack(solver, p, w, c, 1)
    zOptf2 = lex_optimize_knapsack(solver, p, w, c, 2)

    # Compute the ideal point
    zᴵ = [ zOptf1[1] , zOptf2[2] ]

    # Compute one utopian point
    zᵁ = [ zᴵ[1]+1 , zᴵ[2]+1 ]

    # Interactive part
    ρ  = 0.1
    λ = [0.5, 0.5]

    # Compute one nondomainted point in solving the Augmented Weighted Tchebycheff formulation
    z = SolveAugmentedWeightedTchebycheff( solver,  p, w, c,  λ, zᵁ, ρ)

    @printf("  λ=(%6.3f %6.3f)   zND= (%3d %3d) \n",λ[1],λ[2], z[1],z[2])
    return nothing
end


function compute_adaptive_parameters(
            x::Int64, 
            y::Int64
        ) 
    
    if x > y && y ≥ 2
        denom = x*y - y - 3*x + x^2
        w1 = ( x*y - x - y ) / denom
        w2 = ( x*(x - 2)   ) / denom
        ρ  = x / denom
        α  = ( x*y * (x - 1) ) / denom

    elseif y > x && x ≥ 2
        denom = x*y - x - 3*y + y^2
        w1 = ( y * (y - 2) )  / denom
        w2 = (x*y - x - y) / denom
        ρ  = y / denom
        α  = ( x*y * (y - 1) ) / denom   

    elseif x == y && y > 2
        w1 = 0.5
        w2 = 0.5
        ρ  = 1 / (2*(x - 2))
        α  = ( x * (x - 1) ) / ( 2*(x - 2) )
    else
        # Cas dégénéré ou valeurs trop petites
        #w1 = 0.5
        #w2 = 0.5
        #ρ = η * 0.01
        #α = max(x, y)
        error("Cas dégénéré ou valeurs trop petites")
    end

    return w1, w2, ρ
end


# ======================
# Solve subproblem: Augmented Weighted Tchebycheff
# using the Adaptative method
# ======================

function solve_subproblem_knapsack(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64,  
            z1, z2
        )

    rp  = zeros(Int,2);    # reference point
    λ   = zeros(2);        # weights

    x = z2[1]; y = z1[2]   # ideal point
    x = x+1;   y = y+1     # utopian point

    rp[1] = x; rp[2] = y   # setting the reference point with the local utopian point

    print("  Box:  z1=", z1, "  z2=", z2, "  ⇒  RP=($x, $y)   ⟶   ")

    # Compute adaptive parameters
    λ[1], λ[2], ρ = compute_adaptive_parameters(x, y)   

    # Solve the augmented weighted Tchebycheff problem
    z = SolveAugmentedWeightedTchebycheff( solver,  p, w, c,  λ, rp, ρ)

    @printf("λ₁=%6.4f  λ₂=%6.4f   ρ= %8.6f   ⟶   z=(%4d,%4d) \n",  λ[1], λ[2],  ρ, z[1], z[2])

    return z
end



"""
    DachertGorskiKlamroth_knapsack(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64)

    Sequential algorithm for generating the nondominated set based on 
    an adaptive_augmented_tchebycheff. 
    Application to the bi-objective 01 unidimensionnal knapsack problem.  

    In  Kerstin Dächert, Jochen Gorski, Kathrin Klamroth. An augmented weighted 
        Tchebycheff method with adaptively chosen parameters for discrete 
        bicriteria optimization problems. Computers & Operations Research,
        Volume 39, Issue 12, Pages 2929-2943, 2012.
"""
function DachertGorskiKlamroth_knapsack(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64
        )

    # Create two empty sorted lists ND and Temp
    ND   = (Tuple{Int64, Int64})[]
    Temp = (Tuple{Int64, Int64})[]

    # Compute lexicographic optimal solutions z1 and z2 wrt. f1 and f2, respectively
    z1 = lex_optimize_knapsack(solver, p, w, c, 1)
    z2 = lex_optimize_knapsack(solver, p, w, c, 2)

    # Compute the ideal point
    zᴵ = ( z1[1] , z2[2] )

    println( "  z1: ", z1, "   z2: ", z2, "   zᴵ: ", zᴵ)

    if z1 == z2
        # the set of nondominated points is composed of a single point 
        push!(ND, z1)
        return ND
    else
        # Let's play...
        push!(Temp, z1); push!(Temp, z2); #push!(ND, z1)

        while length(Temp) ≥ 2

            # Compute parameters and optionally the new reference point wrt. the first two solutions z1,z2 in Temp
            z1 = Temp[1]; z2 = Temp[2]

            # Solve the resulting subproblem and find solution z
            z  = solve_subproblem_knapsack(solver, p, w, c, z1, z2)

            if z == z1 || z == z2
                # Remove first element from Temp and insert it as last element to ND
                popfirst!(Temp);  push!(ND, z1)
            else
                # Insert z as new second element between z1 and z2 to Temp
                insert!(Temp, 2, z)
            end
        end

        # Add final element of Temp as last element to ND
        push!(ND, Temp[1])

        return ND
    end
end



# =============================================================================
println("-"^80)

solver = GLPK.Optimizer
#solver = HiGHS.Optimizer
#solver = Gurobi.Optimizer
#solver = CPLEX.Optimizer

#p = [3 4 7 9; 5 2 6 4] 
#w = [2, 3, 4, 5]
#c = 8

p = [ 13 10  3 16 12 11  1  9 19 13 ;     # profit 1
       1 10  3 13 12 19 16 13 11  9  ]    # profit 2
w  = [ 4, 4, 3, 5, 5, 3, 2, 3, 5, 4  ]    # weight
c  = 19                                   # capacity


# test SolveAugmentedWeightedTchebycheff ======================================
println("\nAugmented Weighted Tchebycheff:")
rp = [80,80]
ρ  = 0.001
global zPrec=[NaN,NaN]
for step in 0.01:0.01:0.99
    λ =[step,1.0-step] 
    z = SolveAugmentedWeightedTchebycheff( solver,  p, w, c,  λ, rp, ρ)
    if z != zPrec 
        global zPrec = z
        @printf("  λ=(%6.3f %6.3f)   zND= (%3d %3d) \n",λ[1],λ[2], z[1],z[2])
    end
end


# test SteuerChoo =============================================================
println("\nSteuerChoo 1983:")
z = SteuerChoo( solver,  p, w, c)


# test DachertGorskiKlamroth ==================================================
println("\nDachertGorskiKlamroth 2012:")
YN = DachertGorskiKlamroth_knapsack(solver, p, w, c)
println("\n  Set of nondominated points found:")
println("    ", YN)

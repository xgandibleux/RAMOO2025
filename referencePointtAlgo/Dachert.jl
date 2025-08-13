#using JuMP, GLPK 
#using Printf


include("commun.jl")


"""
    compute_adaptive_parameters(
            x::Int64, 
            y::Int64
        )

    Compute the adaptive parameters as described in Table 1 of
    the paper Dächert et al.
"""
function compute_adaptive_parameters(
            x::Int64, 
            y::Int64
        ) 
    
    if x > y && y ≥ 2
        denom = x*y - y - 3*x + x^2
        w1 = ( x*y - x - y ) / denom
        w2 = ( x*(x - 2)   ) / denom
        ρ  = x / denom
        #α  = ( x*y * (x - 1) ) / denom

    elseif y > x && x ≥ 2
        denom = x*y - x - 3*y + y^2
        w1 = ( y * (y - 2) )  / denom
        w2 = (x*y - x - y) / denom
        ρ  = y / denom
        #α  = ( x*y * (y - 1) ) / denom   

    elseif x == y && y > 2
        w1 = 0.5
        w2 = 0.5
        ρ  = 1 / (2*(x - 2))
        #α  = ( x * (x - 1) ) / ( 2*(x - 2) )
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

# ======================

"""
    solve_subproblem_knapsack(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64,  
            z1::Tuple{Int64, Int64}, z2::Tuple{Int64, Int64}
        )

    Solve the Augmented Weighted Tchebycheff subproblem
    using the Adaptative method.
"""
function solve_subproblem_knapsack(  
            solver::DataType, 
            p::Matrix{Int64},  w::Vector{Int64},  c::Int64,  
            z1::Tuple{Int64, Int64}, z2::Tuple{Int64, Int64}
        )

    rp  = zeros(Int,2);    # reference point
    λ   = zeros(2);        # weights

    x = z2[1]; y = z1[2]   # ideal point
    x = x+1;   y = y+1     # utopian point (dominates the ideal point)

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
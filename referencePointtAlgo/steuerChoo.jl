
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

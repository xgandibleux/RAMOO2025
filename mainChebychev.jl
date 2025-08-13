using JuMP, GLPK
using Printf


include("MOCO/knapsack.jl")
include("referencePointtAlgo/commun.jl")


"""
    mainChebychev()

    use SolveAugmentedWeightedTchebycheff for computing YN
"""
function mainChebychev()

    println("-"^80)
    println("\nAugmented Weighted Tchebycheff:")

    solver = GLPK.Optimizer
    p,w,c = initializeMO01UKP()

    rp = [80,80]
    ρ  = 0.001
    zPrec=[NaN,NaN]
    for step in 0.01:0.01:0.99
        λ =[step,1.0-step] 
        z = SolveAugmentedWeightedTchebycheff( solver,  p, w, c,  λ, rp, ρ)
        if z != zPrec 
            zPrec = z
            @printf("  λ=(%6.3f %6.3f)   zND= (%3d %3d) \n",λ[1],λ[2], z[1],z[2])
        end
    end
    
    return nothing
end

mainChebychev()
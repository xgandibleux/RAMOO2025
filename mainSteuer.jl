using JuMP, GLPK
using Printf


include("MOCO/knapsack.jl")
include("referencePointtAlgo/commun.jl")
include("referencePointtAlgo/steuerChoo.jl")


"""
    mainSteuer()

    Use SteuerChoo for computing one nondominated solution interactively
"""
function mainSteuer()

    println("-"^80)
    println("\nSteuerChoo 1983:")

    solver = GLPK.Optimizer
    p,w,c = initializeMO01UKP()

    z = SteuerChoo( solver,  p, w, c)
    return nothing
end

mainSteuer()
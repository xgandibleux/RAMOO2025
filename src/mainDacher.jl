using JuMP, GLPK
using Printf


include("MOCO/knapsack.jl")
include("referencePointtAlgo/commun.jl")
include("referencePointtAlgo/Dachert.jl")


"""
    mainDacher()

    use DachertGorskiKlamroth_knapsack for computing YN
"""
function mainDacher()
    
    println("-"^80)
    println("\nDachertGorskiKlamroth 2012:")

    solver = GLPK.Optimizer
    p,w,c = initializeMO01UKP()

    YN = DachertGorskiKlamroth_knapsack(solver, p, w, c)
    println("\n  Set of nondominated points found:")
    println("    ", YN)

    return nothing
end

mainDacher()
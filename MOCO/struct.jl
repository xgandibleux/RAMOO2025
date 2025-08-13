
"""
Structure pour représenter un point dans l'espace objectif
"""
struct Point2D
    f1::Float64
    f2::Float64
end


Base.:(==)(p1::Point2D, p2::Point2D) = (y1.f1 == y2.f1) && (y1.f2 == y2.f2)
Base.isless(p1::Point2D, p2::Point2D) = y1.f1 < y2.f1 || (y1.f1 == y2.f1 && y1.f2 < y2.f2)

"""
Structure pour représenter un problème d'optimisation bicritère discret
"""
abstract type BicriteriaDiscreteProlem end


"""
Problème de sac à dos biobjectif
"""
mutable struct BiobjectiveKnapsackProblem <: BiobjectiveDiscreteProlem
    p1::Vector{Int}  # coefficients des profits pour objectif 1
    p2::Vector{Int}  # coefficients des profits pour objectif 2
    w::Vector{Int}   # coefficients des poids
    c::Int           # capacite du sac-a-dos
    n::Int           # nombre de variables
    
    function BiobjectiveKnapsackProblem(p1::Vector{Int}, p2::Vector{Int}, w::Vector{Int}, c::Int)
        n = length(p1)
        @assert length(p2) == n && length(a) == n "Toutes les dimensions doivent être égales"
        new(p1, p2, w, c, n)
    end
end

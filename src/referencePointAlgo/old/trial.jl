using LinearAlgebra
using JuMP
using HiGHS

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

"""
Calculer les paramètres optimaux pour la norme Tchebycheff augmentée pondérée
selon les formules dérivées dans le papier (Table 1)
"""
function compute_adaptive_parameters(z1::Point2D, z2::Point2D, η::Float64 = 0.9)
    x = z2.z1
    y = z1.z2
    
    # Cas où x > y >= 2
    if x > y && y >= 2
        denom = x*y - y - 3*x + x^2
        w1 = (x*y - x - y) / denom
        w2 = x*(x - 2) / denom
        ρ = η * x / denom
        α = (x*y*(x - 1)) / denom
        
    # Cas où y > x >= 2
    elseif y > x && x >= 2
        denom = x*y - x - 3*y + y^2
        w1 = y*(y - 2) / denom
        w2 = (x*y - x - y) / denom
        ρ = η * y / denom
        α = (x*y*(y - 1)) / denom
        
    # Cas où x = y > 2
    elseif x == y && x > 2
        w1 = 0.5
        w2 = 0.5
        ρ = η / (2*(x - 2))
        α = (x*(x - 1)) / (2*(x - 2))
        
    else
        # Cas dégénéré ou valeurs trop petites
        w1 = 0.5
        w2 = 0.5
        ρ = η * 0.01
        α = max(x, y)
    end
    
    return w1, w2, ρ, α
end

"""
Résoudre un sous-problème Tchebycheff augmenté pondéré adaptatif
"""
function solve_augmented_weighted_tchebycheff(problem::BicriteriaKnapsackProblem, 
                                             w1::Float64, w2::Float64, ρ::Float64,
                                             reference_point::Point2D,
                                             box_constraints::Union{Nothing, Tuple{Point2D, Point2D}} = nothing)
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    # Variables binaires pour le sac à dos
    @variable(model, x[1:problem.n], Bin)
    
    # Variables pour la linéarisation
    @variable(model, λ >= 0)
    @variable(model, μ1 >= 0)
    @variable(model, μ2 >= 0)
    
    # Objectifs
    @expression(model, z1, sum(problem.c1[i] * x[i] for i in 1:problem.n))
    @expression(model, z2, sum(problem.c2[i] * x[i] for i in 1:problem.n))
    
    # Contrainte du sac à dos
    @constraint(model, sum(problem.a[i] * x[i] for i in 1:problem.n) <= problem.b)
    
    # Contraintes pour la norme Tchebycheff augmentée
    @constraint(model, λ <= w1 * μ1)
    @constraint(model, λ <= w2 * μ2)
    @constraint(model, μ1 <= z1 - reference_point.z1)
    @constraint(model, μ1 <= -(z1 - reference_point.z1))
    @constraint(model, μ2 <= z2 - reference_point.z2)
    @constraint(model, μ2 <= -(z2 - reference_point.z2))
    
    # Contraintes de boîte si spécifiées
    if box_constraints !== nothing
        z1_min, z2_max = box_constraints
        @constraint(model, z1 >= z1_min.z1)
        @constraint(model, z2 <= z2_max.z2)
    end
    
    # Fonction objectif
    @objective(model, Min, λ + ρ * (μ1 + μ2))
    
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        x_sol = value.(x)
        z1_val = sum(problem.c1[i] * x_sol[i] for i in 1:problem.n)
        z2_val = sum(problem.c2[i] * x_sol[i] for i in 1:problem.n)
        return Point2D(z1_val, z2_val), x_sol
    else
        return nothing, nothing
    end
end


"""
Calculer le point idéal global
"""
function compute_ideal_point(problem::BicriteriaKnapsackProblem)
    # Minimiser z1
    model1 = Model(HiGHS.Optimizer)
    set_silent(model1)
    @variable(model1, x[1:problem.n], Bin)
    @constraint(model1, sum(problem.a[i] * x[i] for i in 1:problem.n) <= problem.b)
    @objective(model1, Max, sum(problem.c1[i] * x[i] for i in 1:problem.n))
    optimize!(model1)
    z1_ideal = objective_value(model1)
    
    # Minimiser z2
    model2 = Model(HiGHS.Optimizer)
    set_silent(model2)
    @variable(model2, x[1:problem.n], Bin)
    @constraint(model2, sum(problem.a[i] * x[i] for i in 1:problem.n) <= problem.b)
    @objective(model2, Max, sum(problem.c2[i] * x[i] for i in 1:problem.n))
    optimize!(model2)
    z2_ideal = objective_value(model2)
    
    return Point2D(z1_ideal, z2_ideal)
end

"""
Calculer les solutions lexicographiques minimales
"""
function compute_lexicographic_solutions(problem::BicriteriaKnapsackProblem)
    # Solution lexmin pour f1
    model1 = Model(HiGHS.Optimizer)
    set_silent(model1)
    @variable(model1, x[1:problem.n], Bin)
    @constraint(model1, sum(problem.a[i] * x[i] for i in 1:problem.n) <= problem.b)
    @objective(model1, Min, sum(problem.c1[i] * x[i] for i in 1:problem.n))
    optimize!(model1)
    z1_min = objective_value(model1)
    
    # Parmi les solutions optimales pour f1, minimiser f2
    @constraint(model1, sum(problem.c1[i] * x[i] for i in 1:problem.n) <= z1_min)
    @objective(model1, Min, sum(problem.c2[i] * x[i] for i in 1:problem.n))
    optimize!(model1)
    
    x1_sol = value.(x)
    z1_lex = Point2D(z1_min, objective_value(model1))
    
    # Solution lexmin pour f2
    model2 = Model(HiGHS.Optimizer)
    set_silent(model2)
    @variable(model2, x[1:problem.n], Bin)
    @constraint(model2, sum(problem.a[i] * x[i] for i in 1:problem.n) <= problem.b)
    @objective(model2, Min, sum(problem.c2[i] * x[i] for i in 1:problem.n))
    optimize!(model2)
    z2_min = objective_value(model2)
    
    # Parmi les solutions optimales pour f2, minimiser f1
    @constraint(model2, sum(problem.c2[i] * x[i] for i in 1:problem.n) <= z2_min)
    @objective(model2, Min, sum(problem.c1[i] * x[i] for i in 1:problem.n))
    optimize!(model2)
    
    x2_sol = value.(x)
    z2_lex = Point2D(objective_value(model2), z2_min)
    
    return z1_lex, z2_lex, x1_sol, x2_sol
end

"""
Algorithme séquentiel principal pour générer l'ensemble non-dominé
(Algorithm 1 du papier)
"""
function generate_nondominated_set(problem::BicriteriaKnapsackProblem;
                                  method::Symbol = :adaptive,
                                  use_local_ideal::Bool = false,
                                  η::Float64 = 0.9,
                                  fixed_ρ::Float64 = 1e-3)
    
    println("Calcul des solutions lexicographiques...")
    z1_lex, z2_lex, _, _ = compute_lexicographic_solutions(problem)
    
    if z1_lex == z2_lex
        println("Le point idéal est réalisable")
        return [z1_lex]
    end
    
    # Initialisation
    temp_list = [z1_lex, z2_lex]  # Liste temporaire triée par z1
    sort!(temp_list)
    nd_set = Point2D[]  # Ensemble non-dominé final
    
    ideal_point = compute_ideal_point(problem)
    println("Point idéal: z1=$(ideal_point.z1), z2=$(ideal_point.z2)")
    println("Solutions lexicographiques: $(z1_lex) et $(z2_lex)")
    
    iteration = 0
    
    while length(temp_list) >= 2
        iteration += 1
        println("\n--- Itération $iteration ---")
        
        # Prendre les deux premiers points de temp_list
        z1_current = temp_list[1]
        z2_current = temp_list[2]
        
        println("Exploration entre $(z1_current) et $(z2_current)")
        
        # Choisir le point de référence
        reference_point = if use_local_ideal
            Point2D(min(z1_current.z1, z2_current.z1), min(z1_current.z2, z2_current.z2))
        else
            ideal_point
        end
        
        # Résoudre le sous-problème selon la méthode choisie
        new_point = nothing
        
        if method == :adaptive
            # Méthode adaptative
            w1, w2, ρ, α = compute_adaptive_parameters(z1_current, z2_current, η)
            println("Paramètres adaptatifs: w1=$w1, w2=$w2, ρ=$ρ")
            
            new_point, _ = solve_augmented_weighted_tchebycheff(
                problem, w1, w2, ρ, reference_point,
                use_local_ideal ? (z1_current, z2_current) : nothing
            )
            
        elseif method == :lexicographic
            # Méthode lexicographique
            w1 = z2_current.z2 / (z1_current.z1 + z2_current.z2)
            w2 = z1_current.z1 / (z1_current.z1 + z2_current.z2)
            println("Poids lexicographiques: w1=$w1, w2=$w2")
            
            new_point, _ = solve_lexicographic_tchebycheff(problem, w1, w2, reference_point)
            
        elseif method == :fixed
            # Méthode avec ρ fixe
            w1 = z2_current.z2 / (z1_current.z1 + z2_current.z2)
            w2 = z1_current.z1 / (z1_current.z1 + z2_current.z2)
            println("Poids fixes: w1=$w1, w2=$w2, ρ=$fixed_ρ")
            
            new_point, _ = solve_augmented_weighted_tchebycheff(
                problem, w1, w2, fixed_ρ, reference_point,
                use_local_ideal ? (z1_current, z2_current) : nothing
            )
        end
        
        if new_point === nothing
            println("Aucune solution trouvée")
            popfirst!(temp_list)
            push!(nd_set, z1_current)
        elseif new_point == z1_current || new_point == z2_current
            println("Point existant trouvé: $new_point")
            # Aucun nouveau point dans cette région
            popfirst!(temp_list)
            push!(nd_set, z1_current)
        else
            println("Nouveau point trouvé: $new_point")
            # Insérer le nouveau point dans temp_list
            popfirst!(temp_list)  # Retirer z1_current
            
            # Insérer new_point à sa place ordonnée
            inserted = false
            for i in 1:length(temp_list)
                if new_point.z1 < temp_list[i].z1
                    insert!(temp_list, i, new_point)
                    inserted = true
                    break
                end
            end
            if !inserted
                push!(temp_list, new_point)
            end
            
            push!(nd_set, z1_current)
        end
        
        println("État actuel - temp_list: $(length(temp_list)), nd_set: $(length(nd_set))")
        
        if iteration > 100  # Protection contre les boucles infinies
            println("Limite d'itérations atteinte")
            break
        end
    end
    
    # Ajouter le dernier élément de temp_list à nd_set
    if length(temp_list) == 1
        push!(nd_set, temp_list[1])
    end
    
    sort!(nd_set)
    println("\n=== Résultat final ===")
    println("$(length(nd_set)) points non-dominés trouvés:")
    for (i, point) in enumerate(nd_set)
        println("  $i: $(point)")
    end
    
    return nd_set
end



# Exemple d'utilisation
#function example_usage()
    println("=== Exemple d'utilisation de la méthode Tchebycheff adaptative ===\n")
    
    # Créer un petit problème de test
    problem = BicriteriaKnapsackProblem(
        [13, 10,  3, 16, 12, 11,  1,  9, 19, 13],      # coûts objectif 1
        [1, 10,  3, 13, 12, 19, 16, 13, 11,  9],      # coûts objectif 2
        [4, 4, 3, 5, 5, 3, 2, 3, 5, 4 ],      # profits
        19                     # seuil minimum
    )

    
    println("Problème de sac à dos bicritère avec $(problem.n) items")
    println("Seuil minimum de profit: $(problem.b)")
    
    # Tester différentes méthodes
    methods = [
        (:adaptive, "Méthode Adaptative"),
        (:fixed, "Méthode avec ρ fixe")
    ]
    
    for (method, name) in methods
        println("\n" * "="^50)
        println("$name")
        println("="^50)
        
        time_taken = @elapsed begin
            nd_points = generate_nondominated_set(problem, method=method)
        end
        
        println("Temps de calcul: $(round(time_taken, digits=3)) secondes")
        println("Nombre de points non-dominés: $(length(nd_points))")
    end

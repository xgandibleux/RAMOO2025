using JuMP
using MultiObjectiveAlgorithms
using MathOptInterface
const MOI = MathOptInterface

function analyze_multiobjective_model(model::Model)

    
    # Obtenir le backend MOI du modèle JuMP
    backend = JuMP.backend(model)
    
    #=
    
    # 3. Analyse des objectifs multiples
    println("3. OBJECTIFS MULTIPLES")
    println("-" ^ 40)
    
    # Vérifier si le modèle a des objectifs multiples

    num_objectives = MOI.get(backend, MOI.NumberOfObjectives())
    println("Nombre d'objectifs : $num_objectives")
    
    for obj_idx in 1:num_objectives
        println("\nObjectif $obj_idx:")
        obj_func = MOI.get(backend, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(obj_idx))
        
        println("  Sens: $(MOI.get(backend, MOI.ObjectiveSense(obj_idx)))")
        println("  Terme constant: $(obj_func.constant)")
        println("  Coefficients:")
        
        for term in obj_func.terms
            var_idx = term.variable
            coeff = term.coefficient
            var_name = MOI.get(backend, MOI.VariableName(), var_idx)
            if isempty(var_name)
                var_name = "x$(var_idx.value)"
            end
            println("    $coeff * $var_name")
        end
    end
    println()
    =#
    
    # 4. Analyse des contraintes
    println("4. CONTRAINTES")

    
    # 5. Informations sur les types de variables (entières, binaires, etc.)
    println("5. TYPES DE VARIABLES")
    println("-" ^ 40)
    
    # Variables entières
    if MOI.ConstraintIndex{MOI.VariableIndex, MOI.Integer} in [ci.value for ci in MOI.get(backend, MOI.ListOfConstraints())]
        integer_constraints = MOI.get(backend, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.Integer}())
        if !isempty(integer_constraints)
            println("Variables entières:")
            for constraint in integer_constraints
                var_idx = MOI.get(backend, MOI.ConstraintFunction(), constraint)
                var_name = MOI.get(backend, MOI.VariableName(), var_idx)
                if isempty(var_name)
                    var_name = "x$(var_idx.value)"
                end
                println("  $var_name")
            end
        end
    end
    
    # Variables binaires
    if MOI.ConstraintIndex{MOI.VariableIndex, MOI.ZeroOne} in [ci.value for ci in MOI.get(backend, MOI.ListOfConstraints())]
        binary_constraints = MOI.get(backend, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.ZeroOne}())
        if !isempty(binary_constraints)
            println("Variables binaires:")
            for constraint in binary_constraints
                var_idx = MOI.get(backend, MOI.ConstraintFunction(), constraint)
                var_name = MOI.get(backend, MOI.VariableName(), var_idx)
                if isempty(var_name)
                    var_name = "x$(var_idx.value)"
                end
                println("  $var_name")
            end
        end
    end
    
    println("\n=== FIN DE L'ANALYSE ===")
end

# =============================================================================
function create_multiobjective_model()
    model = Model()
    
    @variable(model, x >= 0)
    @variable(model, y >= 0, Int)
    @variable(model, z, Bin)
    @objective(model, Max, [x + 2*y, 3*x + z])
    @constraint(model, x + 2*y <= 10)
    @constraint(model, 3*x + y >= 5)
    @constraint(model, x + z <= 7)
    
    return model
end




example_model = create_multiobjective_model()
analyze_multiobjective_model(example_model)
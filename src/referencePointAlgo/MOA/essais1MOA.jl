using JuMP
using MathOptInterface
const MOI = MathOptInterface


function afficher_coefficients(model::Model)

    backend = JuMP.backend(model)

    #
    # -------------------------------------------------------------------------
    #
    println("=== VARIABLES ===")
    println("-" ^ 40)    
    num_variables = MOI.get(backend, MOI.NumberOfVariables())
    println("Nombre de variables : $num_variables")

    vars = all_variables(model)
    for (i, v) in enumerate(vars)
        println("[$i] ", name(v))
    end

    #
    # -------------------------------------------------------------------------
    #    
    println("\n=== CONTRAINTES ===")
    println("-" ^ 40)

    # Types de contraintes dans le modèle
    constraint_types = MOI.get(backend, MOI.ListOfConstraintTypesPresent())
    println("Types de contraintes présents : $constraint_types")
    println()    
    
    constraint_count = 1
    for (F, S) in constraint_types
        if F != MOI.VariableIndex  # Ignorer les contraintes de variables (déjà traitées)
            constraints = MOI.get(backend, MOI.ListOfConstraintIndices{F, S}())
            
            for constraint in constraints
                println("Contrainte $constraint_count (Type: $F in $S):")
                
                # Obtenir la fonction de contrainte
                constraint_func = MOI.get(backend, MOI.ConstraintFunction(), constraint)
                constraint_set = MOI.get(backend, MOI.ConstraintSet(), constraint)
                
                # Afficher les coefficients selon le type de fonction
                if isa(constraint_func, MOI.ScalarAffineFunction)
                    println("  Terme constant: $(constraint_func.constant)")
                    println("  Coefficients:")
                    for term in constraint_func.terms
                        var_idx = term.variable
                        coeff = term.coefficient
                        var_name = MOI.get(backend, MOI.VariableName(), var_idx)
                        if isempty(var_name)
                            var_name = "x$(var_idx.value)"
                        end
                        println("    $coeff * $var_name")
                    end
                elseif isa(constraint_func, MOI.VectorAffineFunction)
                    println("  Termes constants: $(constraint_func.constants)")
                    println("  Coefficients (par composante):")
                    for term in constraint_func.terms
                        var_idx = term.scalar_term.variable
                        coeff = term.scalar_term.coefficient
                        output_idx = term.output_index
                        var_name = MOI.get(backend, MOI.VariableName(), var_idx)
                        if isempty(var_name)
                            var_name = "x$(var_idx.value)"
                        end
                        println("    Composante $output_idx: $coeff * $var_name")
                    end
                end
                
                println("  Ensemble de contrainte: $constraint_set")
                println()
                constraint_count += 1
            end
        end
    end

    #
    # -------------------------------------------------------------------------
    #
    println("\n=== OBJECTIFS ===")
    println("-" ^ 40)
    objs = JuMP.objective_function(exemple_modele)
    println("Nombre d'objectifs :", length(objs))

    for (i, obj_expr) in enumerate(objs)
        println("\nObjectif $i :")
        println("  Expression : ", obj_expr)
    end

    vars = all_variables(model)

    for (i, obj_expr) in enumerate(JuMP.objective_function(model))
        println("\nObjectif $i :")
        for v in vars
            coef = JuMP.coefficient(obj_expr, v)
            if coef != 0
                println("  ", name(v), " => ", coef)
            end
        end
    end


end

function create_01UKP()
    model = Model()
    
    @variable(model, x[1:3], Bin)  
    @objective(model, Max, [7*x[1] + 3*x[2], 6*x[2] + 2*x[3]])
    @constraint(model, 2*x[1] + 4*x[2] + x[3] <= 6)

    return model
end

function create_multiobjective_model()
    model = Model()
    @variable(model, x >= 0)
    @variable(model, y >= 0, Int)
    @variable(model, z, Bin)
    @variable(model, 2 <= xv[1:3] <= 10, Int)    
    @objective(model, Max, [x + 2*y + xv[1], 3*x + z])
    @constraint(model, x + 2*y <= 10)
    @constraint(model, 3*x + y >= 5)
    @constraint(model, x + z <= 7)
    return model
end
exemple_modele = create_multiobjective_model()
exemple_modele = create_01UKP()

afficher_coefficients(exemple_modele)

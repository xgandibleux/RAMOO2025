# =============================================================================
function initializeMO01UKP()

    #p = [3 4 7 9; 5 2 6 4] 
    #w = [2, 3, 4, 5]
    #c = 8

    p = [ 13 10  3 16 12 11  1  9 19 13 ;     # profit 1
           1 10  3 13 12 19 16 13 11  9  ]    # profit 2
    w  = [ 4, 4, 3, 5, 5, 3, 2, 3, 5, 4  ]    # weight
    c  = 19                                   # capacity

    return p, w, c
end


# =============================================================================
"""
    generate_MO01UKP(n = 10, o = 2, max_ci = 100, max_wi = 30)

    Generate randomly an instance for the MO-01UKP with c and w uniformly distributed
"""
function generate_MO01UKP(n = 10, o = 2, max_ci = 100, max_wi = 30)

    p = rand(1:max_ci,o,n) # c_i \in [1,max_ci]   # profits
    w = rand(1:max_wi,n) # w_i \in [1,max_wi]     # weight
    c = round(Int64, sum(w)/2)                    # capacity
                
    return p, w, c
end
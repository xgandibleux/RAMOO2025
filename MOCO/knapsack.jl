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

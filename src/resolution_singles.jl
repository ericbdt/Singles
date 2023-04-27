# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX

#include("generation.jl")
include("io_singles.jl")
TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function Singles_cplexSolve(inputFile::String)


    Game = readInputFile(inputFile)
    n = size(Game,1)




    # Create the model
    m = Model(CPLEX.Optimizer)

    @variable(m, M[1:n,1:n], Bin)

    for k in 1:n
        for i in 1:n
            @constraint(m, sum(M[i,j] for j in 1:n if Game[i,j]==k) <= 1)
        end
    end

    for k in 1:n
        for j in 1:n
            @constraint(m, sum(M[i,j] for i in 1:n if Game[i,j]==k) <= 1)
        end
    end    

    ###### Contraintes d'isolement des carrés noirs


    @constraint(m, M[1,2]+M[1,1]>=1)
    @constraint(m, M[2,1]+ M[1,1]>=1)

    @constraint(m, M[1,n-1]+M[1,n]>=1)
    @constraint(m, M[2,n]+ M[1,n]>=1)

    @constraint(m, M[n,2]+ M[n,1]>=1)
    @constraint(m, M[n-1,1]+M[n,1]>=1)

    @constraint(m, M[n,n-1]+ M[n,n]>=1)
    @constraint(m, M[n-1,n]+ M[n,n]>=1)

   for i in 2:n-1
        @constraint(m, M[i-1,1]+ M[i,1]>=1)
        @constraint(m, M[i+1,1]+ M[i,1]>=1)
        @constraint(m, M[i,2]+ M[i,1]>=1)

        @constraint(m, M[i-1,n]+ M[i,n]>=1)
        @constraint(m, M[i+1,n]+ M[i,n]>=1)
        @constraint(m, M[i,n-1]+ M[i,n]>=1)

        @constraint(m, M[1,i-1]+ M[1,i]>=1)
        @constraint(m, M[1,i+1]+ M[1,i]>=1)
        @constraint(m, M[2,i]+ M[1,i]>=1)

        @constraint(m, M[n,i-1]+ M[n,i]>=1)
        @constraint(m, M[n,i+1]+ M[n,i]>=1)
        @constraint(m, M[n-1,i]+ M[n,i]>=1)

        for j in 2:n-1

            @constraint(m, M[i-1,j]+ M[i,j]>=1)
            @constraint(m, M[i+1,j]+ M[i,j]>=1)
            @constraint(m, M[i,j-1]+ M[i,j]>=1)
            @constraint(m, M[i,j+1]+ M[i,j]>=1)

        end
   end

   ##### Contraintes de non isolement des carrés blancs

   @constraint(m, M[1,2]+M[2,1]>= M[1,1])
   @constraint(m, M[n-1,2]+M[n,2]>= M[n,1])
   @constraint(m, M[1,n-1]+M[2,n]>= M[1,n])
   @constraint(m, M[n-1,n]+M[n,n-1]>= M[n,n])
   

    for i in 2:n-1
       @constraint(m, M[i-1,1]+ M[i+1,1] + M[i,2] >= M[i,1])
       @constraint(m, M[i-1,n]+ M[i+1,n] + M[i,n-1] >= M[i,n])
       @constraint(m, M[1,i-1]+ M[1,i+1] + M[2,i] >= M[1,i])
       @constraint(m, M[n,i-1]+ M[n,i+1] + M[n-1,i] >= M[n,i])

        for j in 2:n-1
            @constraint(m, M[i-1,j]+M[i+1,j]+M[i,j-1]+M[i,j+1] >= M[i,j])
        end
    end


    @objective(m,Min,n )




    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(m)


    Ms = JuMP.value.(M)
    Ms = Ms .* Game
    for i in 1:n
        println(Ms[i,:])
    end
        # Return:
    # 1 - true if an optimum is found
    # 2 - the resolution time
    return JuMP.primal_status(m) == JuMP.FEASIBLE_POINT, time() - start
    
end

"""
Heuristically solve an instance
"""
function noircir(i::Int64,j::Int64, M::Matrix{Int64}, cases_vides::Int64, n::Int64, G::Matrix{Int64}) 
    
    new_M = copy(M)
    new_M[i,j] = 0
    new_cv = cases_vides
    new_cv = new_cv - 1

    if i+1 <= n      #la case noire est isolée
        if new_M[i+1,j]==2    
            new_M[i+1,j] = 1
            new_cv = new_cv- 1



            g = G[i+1,j]        ##### On check si on n'a pas engendré deux cases blanches de même valeur sur une même ligne/colonne
            sum = 0
            for k in 1:n 
                if G[i+1,k] == g && new_M[i+1,k]==1
                    sum = sum+1
                end
                if G[k,j] == g && new_M[k,j] == 1
                    sum =sum+1
                end
            end
            if sum > 2
                return true, M, cases_vides
            end
        elseif new_M[i+1,j] == 0
            return true, M, cases_vides #erreur si il y a une noire a cote
        end
    end
    if j+1 <= n 
        if new_M[i,j+1]==2   
            new_M[i,j+1] = 1
            new_cv = new_cv- 1

            g = G[i,j+1]        ##### On check si on n'a pas engendré deux cases blanches de même valeur sur une même ligne/colonne
            sum = 0
            for k in 1:n 
                if G[i,k] == g && new_M[i,k]==1
                    sum = sum+1
                end
                if G[k,j+1] == g && new_M[k,j+1] == 1
                    sum =sum+1
                end
            end
            if sum > 2
                return true, M, cases_vides
            end
        elseif new_M[i,j+1] == 0
            return true, M, cases_vides
        end
    end
    if i-1 >= 1 
        if new_M[i-1,j]==2    
            new_M[i-1,j] = 1
            new_cv = new_cv- 1

            g = G[i-1,j]        ##### On check si on n'a pas engendré deux cases blanches de même valeur sur une même ligne/colonne
            sum = 0
            for k in 1:n 
                if G[i-1,k] == g && new_M[i-1,k]==1
                    sum = sum+1
                end
                if G[k,j] == g && new_M[k,j] == 1
                    sum =sum+1
                end
            end
            if sum > 2
                return true, M, cases_vides
            end
        elseif new_M[i-1,j] == 0
            return true, M, cases_vides
        end
    end
    if j-1 >= 1 
        if new_M[i,j-1]==2   
            new_M[i,j-1] = 1
            new_cv = new_cv- 1

            g = G[i,j-1]        ##### On check si on n'a pas engendré deux cases blanches de même valeur sur une même ligne/colonne
            sum = 0
            for k in 1:n 
                if G[i,k] == g && new_M[i,k]==1
                    sum = sum+1
                end
                if G[k,j-1] == g && new_M[k,j-1] == 1
                    sum =sum+1
                end
            end
            if sum > 2
                return true, M, cases_vides
            end
        elseif new_M[i,j-1] == 0
            return true, M, cases_vides
        end
    end

    for p in 1:n           # on parcourt pour trouver les cases blanches
        for t in 1:n
            if new_M[p,t] == 2
                white = true
                for k in 1:n 
                    if k != t && G[p,t] == G[p,k]
                        if new_M[p,k] > 0   # on verifie que les cases ne sont pas noircies
                            white = false
                        end
                    end 
                    if k != p && G[p,t] == G[k,t]
                        if new_M[k,t] > 0   # on verifie que les cases ne sont pas noircies
                            white = false
                        end
                    end 
                end
                if white 
                    new_M[p,t] = 1
                    new_cv = new_cv - 1
                end
            end 
        end
    end

    for p in 1:n
        for t in 1:n
            if new_M[p,t]==2
                black = false
                for k in 1:n
                    if k != p 
                        if G[p,t] == G[k,t] && new_M[k,t] == 1
                            black = true
                        end
                    end
                    if k != t 
                        if G[p,t] == G[p,k] && new_M[p,k] == 1
                            black = true
                        end
                    end 
                end
                if black 
                    pb, new_M, new_cv = noircir(p,t,new_M,new_cv,n,G)
                    if pb
                        return true, M, cases_vides
                    end
                end
            end
        end
        
    end
    #on teste si les cases blanches ont toutes au moins un voisin blanc 
    #(condition necessaire mais pas suffisante de la convexité)
    for p in 1:n 
        for t in 1:n 
            if new_M[p,t] == 1
                somme = 0 # somme des valeurs voisines
                if p-1 >= 1
                    somme = somme + new_M[p-1,t]
                end
                if t-1 >= 1
                    somme = somme + new_M[p,t-1]
                end
                if p+1 <= n
                    somme = somme + new_M[p+1,t]
                end
                if t+1 <= n
                    somme = somme + new_M[p,t+1]
                end
                if somme == 0
                    return true, M, cases_vides
                end
            end
        end
    end

    return false, new_M, new_cv


end



function heuristicSolve(fichier::String)

    #si on a deja G et n 
    G = readInputFile(fichier)
    n = size(G)[1]
    M = Matrix{Int64}(ones(n,n))
    M = 2 .* M
    cases_vides = n*n

    while cases_vides >0


        #### on parcourt pour trouver les cases blanches initalement détectables
        for i in 1:n           
            for j in 1:n
                if M[i,j] == 2
                    white = true
                    for k in 1:n 
                        if k != j && G[i,j] == G[i,k]
                            if M[i,k] > 0   # on verifie que les cases ne sont pas noircies
                                white = false
                            end
                        end 
                        if k != i && G[i,j] == G[k,j]
                            if M[k,j] > 0   # on verifie que les cases ne sont pas noircies
                                white = false
                            end
                        end 
                    end
                    if white 
                        M[i,j] = 1
                        cases_vides = cases_vides - 1
                    end
                end 
            end
        end



        ##### On cherche maintenant à noircir les cases pouvant l'être

        black = false 
        i=1
        j=1
        while !black && i <= n-1
            while !black && j <= n-1
                if M[i,j] == 2
                    if G[i,j] == G[i,j+1]   ###Notre postulat est qu'il existe toujours au moins deux cases voisines de même nombre
                        if M[i,j+1] >0
                            black = true 
                            pb, M, cases_vides = noircir(i,j,M,cases_vides,n,G)
                            if pb #on teste si cela a créé un probleme, si oui on noircit l'autre
                                pb, M , cases_vides = noircir(i,j+1,M,cases_vides,n,G)
                            end
                        end
                    elseif G[i,j] == G[i+1,j]
                        black = true 
                        pb, M, cases_vides = noircir(i,j,M,cases_vides,n,G)
                        if pb 
                            pb, M, cases_vides= noircir(i+1,j,M,cases_vides,n,G)
                        end
                    end
                end 
                j = j+1
            end
            j = 1
            i = i+1
           
        end
    end

    for i in 1:n
        println((G .* M)[i,:])
    end
end 

"""
Solve all the instances contained in "../data" through CPLEX and heuristics

The results are written in "../res/cplex" and "../res/heuristic"

Remark: If an instance has previously been solved (either by cplex or the heuristic) it will not be solved again
"""
function solveDataSet()

    dataFolder = "../data/"
    resFolder = "../res/"

    # Array which contains the name of the resolution methods
    resolutionMethod = ["cplex"]
    #resolutionMethod = ["cplex", "heuristique"]

    # Array which contains the result folder of each resolution method
    resolutionFolder = resFolder .* resolutionMethod

    # Create each result folder if it does not exist
    for folder in resolutionFolder
        if !isdir(folder)
            mkdir(folder)
        end
    end
            
    global isOptimal = false
    global solveTime = -1

    # For each instance
    # (for each file in folder dataFolder which ends by ".txt")
    for file in filter(x->occursin(".txt", x), readdir(dataFolder))  
        
        println("-- Resolution of ", file)
        readInputFile(dataFolder * file)

        # TODO
        println("In file resolution.jl, in method solveDataSet(), TODO: read value returned by readInputFile()")
        
        # For each resolution method
        for methodId in 1:size(resolutionMethod, 1)
            
            outputFile = resolutionFolder[methodId] * "/" * file

            # If the instance has not already been solved by this method
            if !isfile(outputFile)
                
                fout = open(outputFile, "w")  

                resolutionTime = -1
                isOptimal = false
                
                # If the method is cplex
                if resolutionMethod[methodId] == "cplex"
                    
                    # TODO 
                    println("In file resolution.jl, in method solveDataSet(), TODO: fix cplexSolve() arguments and returned values")
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime = cplexSolve()
                    
                    # If a solution is found, write it
                    if isOptimal
                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write cplex solution in fout") 
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Start a chronometer 
                    startingTime = time()
                    
                    # While the grid is not solved and less than 100 seconds are elapsed
                    while !isOptimal && resolutionTime < 100
                        
                        # TODO 
                        println("In file resolution.jl, in method solveDataSet(), TODO: fix heuristicSolve() arguments and returned values")
                        
                        # Solve it and get the results
                        isOptimal, resolutionTime = heuristicSolve()

                        # Stop the chronometer
                        resolutionTime = time() - startingTime
                        
                    end

                    # Write the solution (if any)
                    if isOptimal

                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write the heuristic solution in fout")
                        
                    end 
                end

                println(fout, "solveTime = ", resolutionTime) 
                println(fout, "isOptimal = ", isOptimal)
                
                # TODO
                println("In file resolution.jl, in method solveDataSet(), TODO: write the solution in fout") 
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end

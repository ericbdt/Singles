# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX

include("generation.jl")
include("io_singles.jl")
TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function Singles_cplexSolve(inputFile::String)


    Game = readInputFile(inputFile)
    n = size(Game,1)




    # Create the model
    m = Model(with_optimizer(CPLEX.Optimizer))

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


    @constraint(m, M[1,2]=1 ,if M[1,1]==0)
    @constraint(m, M[2,1]=1 ,if M[1,1]==0)

    @constraint(m, M[1,n-1]=1 ,if M[1,n]==0)
    @constraint(m, M[2,n]=1 ,if M[1,n]==0)

    @constraint(m, M[n,2]=1 ,if M[n,1]==0)
    @constraint(m, M[n-1,1]=1 ,if M[n,1]==0)

    @constraint(m, M[n,n-1]=1 ,if M[n,n]==0)
    @constraint(m, M[n-1,n]=1 ,if M[n,n]==0)

   for i in 2:n-1
        @constraint(m, M[i-1,1]=1 ,if M[i,1]==0)
        @constraint(m, M[i+1,1]=1 ,if M[i,1]==0)
        @constraint(m, M[i,2]=1 ,if M[i,1]==0)

        @constraint(m, M[i-1,n]=1 ,if M[i,n]==0)
        @constraint(m, M[i+1,n]=1 ,if M[i,n]==0)
        @constraint(m, M[i,n-1]=1 ,if M[i,n]==0)

        @constraint(m, M[1,i-1]=1 ,if M[1,i]==0)
        @constraint(m, M[1,i+1]]=1 ,if M[1,i]==0)
        @constraint(m, M[2,i]=1 ,if M[1,i]==0)

        @constraint(m, M[n,i-1]=1, if M[n,i]==0)
        @constraint(m, M[n,i+1]=1, if M[n,i]==0)
        @constraint(m, M[n-1,i]=1, if M[n,i]==0)

        for j in 2:n-1

            @constraint(m, M[i-1,j]=1, if M[i,j]==0)
            @constraint(m, M[i+1,j]=1, if M[i,j]==0)
            @constraint(m, M[i,j-1]=1, if M[i,j]==0)
            @constraint(m, M[i,j+1]=1, if M[i,j]==0)

        end
   end

   ##### Contraintes de non isolement des carrés blancs

   @constraint(m, M[1,2]+M[2,1]>=1 ,if M[1,1]==1)
   @constraint(m, M[n-1,2]+M[n,2]>=1 ,if M[n,1]==1)
   @constraint(m, M[1,n-1]+M[2,n]>=1 ,if M[1,n]==1)
   @constraint(m, M[n-1,n]+M[n,n-1]>=1 ,if M[n,n]==1)
   

    for i in 2:n-1
       @constraint(m, M[i-1,1]+ M[i+1,1] + M[i,2] >=1 ,if M[i,1]==1)
       @constraint(m, M[i-1,n]+ M[i+1,n] + M[i,n-1] >=1 ,if M[i,n]==1)
       @constraint(m, M[1,i-1]+ M[1,i+1] + M[2,i] >=1 ,if M[1,i]==1)
       @constraint(m, M[n,i-1]+ M[n,i+1] + M[n-1,i] >=1 ,if M[n,i]==1)

        for j in 2:n-1
            @constraint(m, M[i-1,j]+M[i+1,j]+M[i,j-1]+M[i,j+1] >=1 ,if M[i,j]==1)
        end
    end







    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(m)

    # Return:
    # 1 - true if an optimum is found
    # 2 - the resolution time
    return JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start
    
end

"""
Heuristically solve an instance
"""
function noircir(i,j, M, cases_vides) #attention faut mettre le type des variables
    
    new_M = copy(M)
    new_M[i,j] = 0
    new_cv = copy(cases_vides)
    pb = false 

    if i+1 <= n && new_M[i+1,j]==2    #la case noire est isolée
        new_M[i+1,j] = 1
        new_cv = new_cv- 1
    end
    if j+1 <= n && new_M[i,j+1]==2
        new_M[i,j+1] = 1
        new_cv = new_cv- 1
    end
    if i-1 >= 1 && new_M[i-1,j]==2
        new_M[i-1,j] = 1
        new_cv = new_cv- 1
    end
    if j-1 >= 1 && new_M[i,j-1]==2
        new_M[i,j-1] = 1
        new_cv = new_cv- 1
    end

    while new_cv >0 #pas le bon critère, rajouter un critere de si on modif rien on stop
        for p in 1:n
            for t in 1:n
                white = true 
                for k in 1:n 
                    if k != t && G[p,t] == G[p,k] #probleme : les valeurs de G ne changent pas, il faudrait les mettre à zero quand M ij passe à 0
                        if new_M[p,t] 
                        white = false
                    end 
                    if k != p && G[p,t]==G[k,t]
                        white = false
                    end 
                end
                if white 
                    new_M[p,t] = 1
                    new_cv = new_cv- 1
                end 
            end
        end

        #parcourir la grille
        #si une case contient un chiffre qui est déjà dans la ligne/colonne
        #avec une valeur 1 dans M 
        #alors cette case est noire 
        #tester les contraintes (tester si il y a une noire a cote)
        if pb
            return pb, [], 1 # la matrice et le nb sont inutiles
        else
            #on la noircit 
        end

        
    end



end



function heuristicSolve()

    #si on a deja G et n 

    M = Matrix{Int64}(ones(n,n))
    M = 2.*M
    cases_vides = n*n

    while cases_vides >0
        for i in 1:n
            for j in 1:n
                white = true
                for k in 1:n 
                    if k != j && G[i,j] == G[i,k]
                        white = false
                    end 
                    if k != i && G[i,j] == G[k,j]
                        white = false
                    end 
                end
                if white 
                    M[i,j] = 1
                    cases_vides = cases_vides -1
                end 
            end
        end

        black = false 
        i=1
        j=1
        while not black && i <= n-1
            while not black && j <= n-1
                if G[i,j] == G[i,j+1]
                    black = true 
                    cases_vides=cases_vides-1
                    pb, new_M, new_cv = noircir(i,j,M)
                    if pb 
                        bool, M , cases_vides= noircir(i,j+1,M)
                    else
                        M = new_M
                        cases_vides = new_cv
                    end
                elseif G[i,j] == G[i+1,j]
                    black = true 
                    cases_vides=cases_vides-1
                    M[i,j] = 0
                    pb, new_M, new_cv = noircir(i,j,M)
                    if pb 
                        M[i,j] = 1
                        M[i+1,j] = 0
                        bool, M, cases_vides= noircir(i+1,j,M)
                    else
                        M = new_M
                        cases_vides = new_cv
                    end

                end 
                j = j+2
            end
            j = 1
            i = i+2
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

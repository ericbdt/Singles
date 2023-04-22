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
function heuristicSolve()

    # TODO
    println("In file resolution.jl, in method heuristicSolve(), TODO: fix input and output, define the model")
    
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

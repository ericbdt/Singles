# This file contains methods to generate a data set of instances (i.e., sudoku grids)
include("io_singles.jl")


function generateInstance(n::Int64)

    G=Matrix{Int64}(ones(n,n))

    nb_black = rand(0:(div(n+1,2)))
    top_black=[]
    for i in 1:n

        #determination des cases noires
        black = []
        sides = []
        for j in 1:nb_black
            k = rand(1:n)
            while k in top_black || k in sides || k in black
                println("1")
                k = rand(1:n)
            end
            println("sorti")
            G[i,k] = 0 
            push!(black,k)
            if k-1 >=1
                push!(sides,k-1)
            end
            if k+1<=n
                push!(sides,k+1)
            end
        end
        
        #valeurs des cases blanches
        for j in 1:n 
            if G[i,j] == 1
                k = rand(1:n)
                while k in G[i,1:(j-1)] || k in G[1:(i-1),j]
                    k = rand(1:n)    
                    println("2")
                end
                G[i,j]=k
            end
        end

        nb_black = rand(0:(div(n+1,2) - nb_black))
        top_black = black

        
    end

    #ajout des valeurs des cases noires

    for i in 1:n 
        for j in 1:n 
            if G[i,j] == 0
                choices = [G[i,:];G[:,j]]
                k = rand(choices)
                while k == 0
                    k = rand(choices)
                    println("3")
                end
                G[i,j] = k
            end
        end
    end


    for i in 1:n 
        println(G[i,:])
    end


    fichier = open("data/singles_taille_"*string(n)*".txt","w")
    line=""
    for k in 1:n
        line=""
        for l in 1:n-1
            line = line*string(G[k,l])
            line = line*","
        end
        line = line*string(G[k,n])
        line = line*"\n"
        write(fichier, line)
    end
    close(fichier)




end 

"""
Generate all the instances

Remark: a grid is generated only if the corresponding output file does not already exist
"""
function generateDataSet()

    # TODO
    println("In file generation.jl, in method generateDataSet(), TODO: generate an instance")
    
end




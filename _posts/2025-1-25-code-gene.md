---
title: genetic algorithm
categories: [code]
comments: true
---


```julia
using Random
randperm(5)

```




    5-element Vector{Int64}:
     5
     1
     2
     3
     4




```julia
struct solution 
    best::Float64
    bestIndividual::Vector{Float64}
    convergence::Vector{Float64}
end
```


```julia
function genetic_algorithm(f::Function, population::Vector{Vector{Float64}}, S, C, M 
                                                    ;  max_gen=10, p_c=0.4, p_m=0.01)
    best_current_convergence = zeros(max_gen)
    convergence  = zeros(max_gen)
    best_so_far_gene= nothing
    best_so_far_fitness=0
            
    for k in 1:max_gen
        fitness = f.(population)
        parents::Vector{Int} = select(S, fitness)
        population = population[parents]
        
        new_population::Vector{Vector{Float64}}=[]
        while (population != [] )
            if (rand()<p_c && length(population)>1) 
                parent1= pop!(population)
                parent2= pop!(population)
                children::Vector{Vector{Float64}} = crossover(C, parent1,parent2)
                push!(new_population, children[1])
                push!(new_population, children[2])
            
            elseif rand()<p_m
                parent1=pop!(population)
                mutant = mutate(M, parent1)
                push!(new_population, mutant)
            else
                child=pop!(population)
                push!(new_population, child)
            end
        end
        population=deepcopy(new_population)
        
        
        best_idx= argmax(f.(population))
        best_current_fitness = f(population[best_idx])
        if (best_so_far_fitness < best_current_fitness)
            best_so_far_fitness=best_current_fitness
            best_so_far_gene=population[best_idx]
        end
        best_current_convergence[k]= best_current_fitness
        convergence[k] = best_so_far_fitness
    
    end
    s = solution(best_so_far_fitness, best_so_far_gene, convergence)
    return s, best_current_convergence
end
```




    genetic_algorithm (generic function with 1 method)




```julia
abstract type SelectionMethod end

###
struct TruncationSelection <: SelectionMethod
    k
end

function select(t::TruncationSelection, y)
    parents_index= sortperm(y, rev=true)
    return [parents_index[rand(1:t.k)] for _ in y] 
end



###
struct TournamentSelection <: SelectionMethod
    k
end

function select(t::TournamentSelection, y)
    pool=[]
    
    n=length(y)
    for _ in 1:n
        parents_index=randperm(n) 
        push!(pool, parents_index[argmax(y[parents_index[1:t.k]])])
    end
    
    return [rand(pool) for _ in y]
end



###
struct StochasticUniversalSampling <: SelectionMethod
end

function select(t::StochasticUniversalSampling, y)
    pool=[]
    n=length(y)
    y_sum=sum(y)
    
    roulette=[]; v=0
    for i in 1:n
        v += y[i]/y_sum
        push!(roulette, v)
    end
    
    region=0
    ptr=rand()/n
    for i in 1:n
        region+= y[i] / y_sum
            
        while(ptr<region)
            parent_index=argmax(roulette .> ptr)
            push!(pool,parent_index)
            ptr+=1/n
        end
    end
    
   return pool         
end
```




    select (generic function with 3 methods)




```julia
abstract type CrossoverMethod end

###
struct SinglePointCrossover <: CrossoverMethod 
end

function crossover(::SinglePointCrossover, parent1, parent2)
    idx=rand(1:length(parent1))
    child1= vcat(parent1[1:idx],parent2[idx+1:end]) 
    child2= vcat(parent2[1:idx],parent1[idx+1:end]) 
    return [child1, child2]
end



###
struct TwoPointCrossover <: CrossoverMethod
end

function crossover(::TwoPointCrossover, parent1, parent2)
    n= length(parent1)
    i,j = rand(1:n, 2)
    if i>j
        (i,j)=(j,i)
    end 
    
    child1=vcat(parent1[1:i],parent2[i+1:j],parent1[j+1:end]) 
    child2=vcat(parent2[1:i],parent1[i+1:j],parent2[j+1:end]) 
    return [child1, child2]
end



###
struct UniformCrossover <: CrossoverMethod
end

function crossover(::UniformCrossover, parent1, parent2)
    child1 = deepcopy(parent1)
    child2 = deepcopy(parent2)
    for i in 1:length(parent1)
        if rand()<0.5
            child1[i]=parent2[i]
            child2[i]=parent1[i]
        end
    end
    
    return [child1, child2]
end

```




    crossover (generic function with 3 methods)




```julia
abstract type MutationMethod end

struct GaussianMutation <: MutationMethod
    σ
end

function mutate(M::GaussianMutation, child)
    return child + randn(length(child)) * M.σ
end
```




    mutate (generic function with 1 method)




```julia
function f(x)
    500-sum(x.^2)
end

ub=10
lb=-10
m=100
population= [rand(Float64,(4))*(ub-lb).+lb for _ in 1:m]

max_gen=1000

S=StochasticUniversalSampling()#TournamentSelection(10)#TruncationSelection(10)#
C=SinglePointCrossover()#UniformCrossover()#TwoPointCrossover()#
M=GaussianMutation(0.1)


s,curve =genetic_algorithm(f,population,S,C,M; max_gen=max_gen)
#@show s
println("best score: ",s.best ,"\n best solution:", s.bestIndividual)

using PyPlot

plot(1:max_gen, 500 .- s.convergence)
plot(1:max_gen, 500 .- curve)
```

    best score: 499.70089949427177
     best solution:[0.3320951436064817, 0.3488988231049739, 0.25884738183063094, 0.008998081830054502]
    


    
![png](/_posts//picture/output_6_1.png)
    





    1-element Vector{PyCall.PyObject}:
     PyObject <matplotlib.lines.Line2D object at 0x0000028B15D1EF00>




```julia
population

```




    5-element Vector{Vector{Float64}}:
     [9.339967677238697, 6.001535739932116, -7.734812969555902, 4.952256486998143]
     [-3.1744404779188606, -3.8244228475920927, 4.1025917149242765, -3.497103024763666]
     [5.549237067904224, 9.016052339010944, -6.883996913084005, -7.281847206246144]
     [7.893339058972796, -1.3767624586748877, 1.8732078544175081, 1.195583133582419]
     [-9.41375332010027, -2.8790902703435606, 8.253411022881515, 9.97455100588843]




```julia

```

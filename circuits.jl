using Random

function stabilizerbell(N::Integer)
    N % 2 == 0 || error("cannot make Bell pairs out of an odd number of sites.")
    circuit = AbstractSymbolicOperator[]
    for i in 1:(N÷2)
        push!(circuit, sHadamard(2i-1))
        push!(circuit, sCNOT(2i-1, 2i))
    end
    return circuit
end

function randomlayerspecs(N::Integer, nlayers::Integer; seed::Union{Nothing,Integer}=nothing, psingle::Float64=0.6, pcnot::Float64=0.5)
    N ≥ 1 || error("the number of qubits must be at least 1.")
    nlayers ≥ 1 || error("the number of layers must be at least 1.")
    0.0 ≤ psingle ≤ 1.0 || error("psingle must be between 0 and 1.")
    0.0 ≤ pcnot ≤ 1.0 || error("pcnot must be between 0 and 1.")

    rng = isnothing(seed) ? Random.default_rng() : Random.MersenneTwister(seed)
    specs = Tuple{Symbol,Int,Int}[]

    for _ in 1:nlayers
        for i in 1:N
            if rand(rng) < psingle
                gate = rand(rng, (:(H), :(S)))
                push!(specs, (gate, i, 0))
            end
      end

        if N ≥ 2
            perm = randperm(rng, N)
            for k in 1:2:(N - 1)
                if rand(rng) < pcnot
                    a = perm[k]
                    b = perm[k + 1]
                    if rand(rng) < 0.5
                        push!(specs, (:(CNOT), a, b))
                    else
                        push!(specs, (:(CNOT), b, a))
                    end
                end
          end
      end
    end

    return specs
end

function stabilizerrandomlayers(N::Integer, nlayers::Integer; seed::Union{Nothing,Integer}=nothing, psingle::Float64=0.6, pcnot::Float64=0.5)
  specs = randomlayerspecs(N, nlayers; seed=seed, psingle=psingle, pcnot=pcnot)
  return stabilizerrandomlayers(specs)
end

function stabilizerrandomlayers(specs)
  circuit = AbstractSymbolicOperator[]
  for (gate, i, j) in specs
        if gate == :H
            push!(circuit, sHadamard(i))
        elseif gate == :S
            push!(circuit, sPhase(i))
        elseif gate == :CNOT
            push!(circuit, sCNOT(i, j))
        end
  end

  return circuit
end

function mpsrandomlayers(N::Integer, nlayers::Integer, sites; seed::Union{Nothing,Integer}=nothing, psingle::Float64=0.6, pcnot::Float64=0.5)
  N == length(sites) || error("N must match the number of provided sites.")
  specs = randomlayerspecs(N, nlayers; seed=seed, psingle=psingle, pcnot=pcnot)
  return mpsrandomlayers(specs, sites)
end

function mpsrandomlayers(specs, sites)
  circuit = ITensor[]
  for (gate, i, j) in specs
        if gate == :H
            push!(circuit, op("H", sites[i]))
        elseif gate == :S
            push!(circuit, op("S", sites[i]))
        elseif gate == :CNOT
            push!(circuit, op("CNOT", sites[i], sites[j]))
        end
  end

  return circuit
end

function mpsbell(N::Integer, sites)
    N % 2 == 0 || error("cannot make Bell pairs out of an odd number of sites.")
    circuit = ITensor[]
    for i in 1:(N÷2)
        # Each push! adds an operator at the end of the circuit
        push!(circuit, op("H", sites[2i-1]))
        push!(circuit, op("CNOT", sites[2i-1], sites[2i]))
    end
    return circuit
end

function stabilizergeneric1(N::Integer)
    N  ≥ 8 || error("the chosen generic circuit must have at least 8 qubits.") 
    circuit = AbstractSymbolicOperator[]

    for i in (1,2,3,4,5,6,7,8)
        push!(circuit, sHadamard(i))
    end
    for (i,j) in zip((1,4,5), (2,6,7))
      push!(circuit, sCNOT(i, j))
    end
    for i in (3,5,6)
      push!(circuit, sHadamard(i))
    end
    for (i,j) in zip((1,4), (6,7))
      push!(circuit, sCNOT(i, j))
    end
    for i in (2,3,5,8)
        push!(circuit, sHadamard(i))
    end
    for (i,j) in zip((1,2,4,5,7), (2,3,6,7,8))
      push!(circuit, sCNOT(i, j))
    end
    for i in (1,5,6,7)
        push!(circuit, sHadamard(i))
    end
    for (i,j) in zip((1,2), (4,8))
      push!(circuit, sCNOT(i, j))
    end
    return circuit
end

function mpsgeneric1(N::Integer, sites)
    N  ≥ 8 || error("the chosen generic circuit must have at least 8 qubits.") 
    circuit = ITensor[]

    for i in (1,2,3,4,5,6,7,8)
        push!(circuit, op("H", sites[i]))
    end
    for (i,j) in zip((1,4,5), (2,6,7))
      push!(circuit, op("CNOT", sites[i], sites[j]))
    end
    for i in (3,5,6)
      push!(circuit, op("H", sites[i]))
    end
    for (i,j) in zip((1,4), (6,7))
      push!(circuit, op("CNOT", sites[i], sites[j]))
    end
    for i in (2,3,5,8)
        push!(circuit, op("H", sites[i]))
    end
    for (i,j) in zip((1,2,4,5,7), (2,3,6,7,8))
      push!(circuit, op("CNOT", sites[i], sites[j]))
    end
    for i in (1,5,6,7)
        push!(circuit, op("H", sites[i]))
    end
    for (i,j) in zip((1,2), (4,8))
      push!(circuit, op("CNOT", sites[i], sites[j]))
    end
    return circuit
end


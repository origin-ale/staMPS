using ITensors, ITensorMPS
using QuantumClifford, Quantikz
using Printf

function stabevolve(N::Integer)
    println("""Running stabilizer simulation on $N qubits starting at |0⟩|0⟩…|0⟩.
        Creating 4 Bell pairs.
        """)

    state = one(Destabilizer, N) # Diagonal Z stabilizer tableau, i.e. start from |0⟩|0⟩…|0⟩

    printstyled("=== Initial tableau ===\n"; color= :blue)
    println(state)

    circuit = AbstractSymbolicOperator[]

    push!(circuit, sHadamard(1))
    push!(circuit, sHadamard(3))
    push!(circuit, sHadamard(5))
    push!(circuit, sHadamard(7))

    push!(circuit, sCNOT(1,2))
    push!(circuit, sCNOT(3,4))
    push!(circuit, sCNOT(5,6))
    push!(circuit, sCNOT(7,8))

    stats = @timed mctrajectory!(state, circuit) # Apply the circuit to the state in-place (one MC iter, but this is deterministic)

    @printf "\nEvolution done in %.3f s.\n" stats.time
    printstyled("\n=== Final tableau ===\n"; color= :blue)
    println(state)
end

function mpsevolve(N::Integer)
    println("""Running MPS simulation on $N qubits starting at |0⟩|0⟩…|0⟩.
        Creating 4 Bell pairs.
        """)

    sites = siteinds("Qubit", N) # Array of one physical index per site
    states = ["↑" for n in 1:N] # Array of single-qubit states
    ψ = MPS(sites, states) # MPS with site indices, in product state |0⟩|0⟩…|0⟩

    circuit = ITensor[]

    # Each push! adds an operator at the end of the circuit
    push!(circuit, op("H", sites[1]))
    push!(circuit, op("H", sites[3]))
    push!(circuit, op("H", sites[5]))
    push!(circuit, op("H", sites[7]))

    push!(circuit, op("CNOT", sites[1], sites[2]))
    push!(circuit, op("CNOT", sites[3], sites[4]))
    push!(circuit, op("CNOT", sites[5], sites[6]))
    push!(circuit, op("CNOT", sites[7], sites[8]))

    ψ, stats... = @timed apply(circuit, ψ)

    @printf "Evolution done in %.3f s.\n" stats.time
end

N=8
stabevolve(N)

println("\n")

mpsevolve(N)

println("\n")
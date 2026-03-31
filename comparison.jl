using ITensors, ITensorMPS
using QuantumClifford, Quantikz
using Printf

function stabilizerbell(N::Integer)
    N % 2 == 0 || error("cannot make Bell pairs out of an odd number of sites.")
    circuit = AbstractSymbolicOperator[]
    for i in 1:(N÷2)
        push!(circuit, sHadamard(2i-1))
        push!(circuit, sCNOT(2i-1, 2i))
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

function stabevolve(N::Integer)
    println("""Running stabilizer simulation on $N qubits starting at |0⟩|0⟩…|0⟩.
        Creating $(N÷2) Bell pairs.
        """)

    state = one(Destabilizer, N) # Diagonal Z stabilizer tableau, i.e. start from |0⟩|0⟩…|0⟩
    circuit = stabilizerbell(N)

    stats = @timed mctrajectory!(state, circuit) # Apply the circuit to the state in-place (one MC iter, but this is deterministic)

    @printf "Evolution done in %.3f s.\n" stats.time
    return state
end

function mpsevolve(N::Integer)
    println("""Running MPS simulation on $N qubits starting at |0⟩|0⟩…|0⟩.
        Creating $(N÷2) Bell pairs.
        """)

    sites = siteinds("Qubit", N) # Array of one physical index per site
    states = ["↑" for n in 1:N] # Array of single-qubit states

    ψ = MPS(sites, states) # MPS with site indices, in product state |0⟩|0⟩…|0⟩
    circuit = mpsbell(N, sites)

    ψ, stats... = @timed apply(circuit, ψ)

    @printf "Evolution done in %.3f s.\n" stats.time
    return ψ
end

function convertstabmpo(stabstr::AbstractString, sites)
    clean = replace(strip(stabstr), " " => "")
    isempty(clean) && error("Encountered empty stabilizer string")

    phasechar = clean[1]
    paulis = clean[2:end]
    length(paulis) == length(sites) || error("Stabilizer length does not match MPS size")

    coeff = phasechar == '+' ? 1.0 : phasechar == '-' ? -1.0 : error("Unsupported stabilizer phase in $stabstr")

    terms = Any[coeff] # First element in opsum terms is the coefficient (sign)
    for (j, p) in enumerate(paulis)
        opname = if p == 'X'
            "X"
        elseif p == 'Y'
            "Y"
        elseif p == 'Z'
            "Z"
        elseif p == '_' || p == 'I'
            "Id"
        else
            error("Unsupported Pauli symbol '$p' in $stabstr")
        end
        
        push!(terms, opname) # Then each non-identity Pauli gets added…
        push!(terms, j) # …along with the site it acts on
    end

    os = OpSum()
    add!(os, terms...)
    return MPO(os, sites)
end

function checkstamps(stabψ::Destabilizer, mpsψ::MPS; atol::Float64=1e-10)
    sites = siteinds(mpsψ)
    stabilizerstrings = [string(s) for s in stabilizerview(stabψ)]

    ψref = deepcopy(mpsψ)
    normalize!(ψref)

    println("Checking compatibility…")
    all_ok = true
    for (k, sstr) in enumerate(stabilizerstrings)
        mpo = convertstabmpo(sstr, sites)
        ψapplied = apply(mpo, ψref)
        normalize!(ψapplied)

        ov = inner(ψref, ψapplied)
        ok = abs2(ov - 1) ≤ atol
        all_ok &= ok

        status = ok ? "PASS" : "FAIL"
        @printf "Stabilizer %d (%s): |⟨ψ|Sψ⟩|² = %.12f [%s]\n" k sstr ov status
    end

    if all_ok
        printstyled("\nFinal states are equal.", color= :green)
      else
        printstyled("\nWARNING: at least one stabilizer check failed. Final states are not equal, check code and rerun!", color= :red)
    end

    return all_ok
end


N = 6

stabψ = stabevolve(N)
println("\n")

mpsψ = mpsevolve(N)
println("\n")

checkstamps(stabψ, mpsψ)
println("\n")
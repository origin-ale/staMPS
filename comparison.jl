include("circuits.jl")

using ITensors, ITensorMPS
using QuantumClifford, Quantikz
using Printf
using ProgressMeter

function stabevolve(N, circuitspecs)
    printstyled("--- Stabilizer simulation ---\n"; color = :cyan)

    state = one(Destabilizer, N) # Diagonal Z stabilizer tableau, i.e. start from |0⟩|0⟩…|0⟩
    circuit = stabilizerrandomlayers(circuitspecs)

    print("Evolving… ")
    stats = @timed mctrajectory!(state, circuit) # Apply the circuit to the state in-place (one MC iter, but this is deterministic)

    @printf "done in %.3f s.\n" stats.time
    return state
end

function mpsevolve(N, circuitspecs)
    printstyled("--- MPS simulation ---\n"; color = :cyan)

    sites = siteinds("Qubit", N) # Array of one physical index per site
    states = ["↑" for n in 1:N] # Array of single-qubit states

    ψ = MPS(sites, states) # MPS with site indices, in product state |0⟩|0⟩…|0⟩
    circuit = mpsrandomlayers(circuitspecs, sites)

    @showprogress desc = "Evolving…" color = :normal for gate in circuit
        ψ = apply(gate, ψ)
    end

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

    printstyled("--- Compatibility check ---\n"; color = :cyan)
    println("Computing action of stabilizers on final MPS…\n")
    all_ok = true
    for (k, sstr) in enumerate(stabilizerstrings)
        mpo = convertstabmpo(sstr, sites)
        ψapplied = apply(mpo, ψref)
        normalize!(ψapplied)

        ov = abs2(inner(ψref, ψapplied))
        ok = abs2(ov - 1) ≤ atol^2
        all_ok &= ok

        if length(stabilizerstrings) ≤ 16 # Avoid flooding terminal for many-body systems
          status = ok ? "PASS" : "FAIL"
          @printf "Stabilizer %d (%s): |⟨ψ|Sψ⟩|² = %.12f [%s]\n" k sstr ov status
        end
    end

    if all_ok
        printstyled("\nFinal states are equal.", color= :green)
      else
        printstyled("\nWARNING: at least one stabilizer check failed. Final states are not equal, check code and rerun!", color= :red)
    end

    return all_ok
end


N = 14
layers = 50
seed = 42

printstyled("===== COMPARING MPS AND STABILIZER FORMALISM =====\n\n", color = :cyan)
circuitspecs = randomlayerspecs(N, layers; seed=seed)

stabψ = stabevolve(N, circuitspecs)
println("\n")

mpsψ = mpsevolve(N, circuitspecs)
println("\n")

checkstamps(stabψ, mpsψ)
println("\n")
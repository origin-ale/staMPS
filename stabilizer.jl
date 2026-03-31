using QuantumClifford, Quantikz

N = 8 # Number of qubits

println("""Running stabilizer simulation on $N qubits starting at |0⟩|0⟩…|0⟩.
        Random all-to-all Clifford circuit.
        """)

state = one(Destabilizer, N) # Diagonal Z stabilizer tableau, i.e. start from |0⟩|0⟩…|0⟩

printstyled("=== Initial state ===\n"; color= :blue)
println(state)

circuit = random_all_to_all_clifford_circuit(N, 3N) # 3N random 2-qubit cliffs acting on the N qubits

mctrajectory!(state, circuit) # Apply the circuit to the state in-place (one MC iter, but this is deterministic)

printstyled("\n=== Final state ===\n"; color= :blue)
println(state)

println('\n')
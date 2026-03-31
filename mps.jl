using ITensorMPS
using ITensors

N = 8 # Number of qubits

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

ψ = apply(circuit, ψ)
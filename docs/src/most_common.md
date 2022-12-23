# [CNOT-like Purification Gates](@id cnotlike-gates)

```@meta
DocTestSetup = quote
    using QuantumClifford, BPGates
end
```

As a reminder, we have the following basis states (written in binary, in decimal, and as list of stabilizers):
```
00 = 0 = +XX
         +ZZ

01 = 1 = +XX
         -ZZ

10 = 2 = -XX
         +ZZ

11 = 3 = -XX
         -ZZ
```

If we consider purification circuits that preserve the `00` state, the majority of useful gates in purification circuits can be written in the form:

```@raw html
<img style="width:50%" src="../a_good_bp_gate.png">
```

Above the entangled pairs are qubit 1&3 and qubits 2&4.

The main property of this type of circuits is that they always map the "good" state to itself, while they permute the rest of the Bell basis states.

The $h_a\otimes h_b$ (and $f_a\otimes f_b$) can be one of the six permutations of `01, 10, 11`, implemented by single-qubit Clifford gates. There are 6 possible such operations.
The hard-coded CNOT gate provides the necessary entangling.

These gates are implemented as [`CNOTPerm`](@ref).

The 6 permutations (i.e. the gates $h_a\otimes h_b$ and $f_a\otimes f_b$) are:

```@setup 1
using BPGates # hide
```

```@repl 1
BPGates.good_perm_tuple
```

The gates represented as Clifford operations are:

TODO finish this listing
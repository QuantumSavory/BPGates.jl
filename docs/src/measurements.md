# [Coincidence Measurements](@id measurements)

```@meta
DocTestSetup = quote
    using QuantumClifford, BPGates
end
```

Coincidence measurements are implemented by [`BellMeasure`](@ref). There are three possible types of measurements, one for each basis `X`, `Y`, and `Z`.

`BellMeasure` is supported by `QuantumClifford`'s `apply!`, but more importantly, it is supported by `mctrajectories`, which provides for an easy way to run multiple Monte Carlo simulation trajectories. A successful coincidence measurement lets the simulations continue. A failed coincidence ends the simulations with a reported error.

```@setup 1
using BPGates
using QuantumClifford
```

### Measurement results depending on basis

As seen in the code example below, the meaning of the measurements is as follows

|Call|Basis & Check|States passed|
|---|---|---|
|`BellMeasure(1,_)`|X coincidence|`00`&`01`|
|`BellMeasure(2,_)`|Y anticoincidence|`00`&`11`|
|`BellMeasure(3,_)`|Z coincidence|`00`&`10`|

For reference, here are the states in different representations

| `BPGates` notation| Stabilizer tableaux | Kets | in X basis | in Y basis |
|:---|:---|:---|:---|:---|
|`00`|`+XX +ZZ`|`∣00⟩+∣11⟩`|`∣++⟩+∣--⟩`|`∣i₊i₋⟩+∣i₋i₊⟩`|
|`01`|`+XX -ZZ`|`∣01⟩+∣10⟩`|`∣++⟩-∣--⟩`|`∣i₊i₊⟩-∣i₋i₋⟩`|
|`10`|`-XX +ZZ`|`∣00⟩-∣11⟩`|`∣+-⟩+∣-+⟩`|`∣i₊i₊⟩+∣i₋i₋⟩`|
|`11`|`-XX -ZZ`|`∣01⟩-∣10⟩`|`∣+-⟩-∣-+⟩`|`∣i₊i₋⟩-∣i₋i₊⟩`|

```@repl 1
all_states = [BellState([0,0]), BellState([1,0]), BellState([0,1]), BellState([1,1])];
filter_true(meas) = [state
    for state in all_states
    if bellmeasure!(copy(state), meas)[2]];
filter_true(BellMeasure(1,1))
filter_true(BellMeasure(2,1))
filter_true(BellMeasure(3,1))
```
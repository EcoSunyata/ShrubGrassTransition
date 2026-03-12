# ShrubGrassTransition

This repository contains the simulation code accompanying the manuscript:
"Phase transition in a dryland shrub–grass ecosystem: Emergence through local trade-offs".

The model explores how simple, local resource-foraging trade-offs between shrubs (extensive strategy) and grasses (intensive strategy) can give rise to complex, system-level phase transitions, including pattern formation, abrupt shifts, boundary effects, and the emergence of fairy circles.

## Model overview

Plants are modeled as self-replicating units that interact via a spatially explicit kernel. Each species (shrub or grass) generates an interaction field based on its kernel, which combines both competitive (negative) and facilitative (positive) effects.

The key trade-off is between:
Shrubs: Extensive root systems, wide negative influence (low intensity), higher survival threshold.
Grasses: Intensive resource use, strong but narrow negative influence, lower survival threshold.

The system evolves based on local interaction fields, leading to emergent phenomena such as:
1. Turing-like spatial patterns
2. Abrupt state transitions
3. Boundary effects and collapse
4. Fairy circle formation

## Notebook Descriptions

### 1. `ShrubGrassTransition.ipynb`

- Implements the kernel-based cellular automaton.
    
- Simulates competition between shrubs and grasses across a resource gradient (`s`).
    
- Produces density dynamics and final spatial patterns.
    
- Includes sensitivity analysis across resource levels.
    

### 2. `BoundaryEffect.ipynb`

- Models the interface between vegetation and desert.
    
- Introduces disturbance and enhanced positive feedback.
    
- Demonstrates how patch expansion and sudden collapse emerge under increasing disturbance.
    

### 3. `FariyCircle.ipynb`

- Simulates fairy circle formation through "kernel contraction" (reduced foraging range).
    
- Shows how shrubs transition into hollow rings under competition from grasses.
    
- Demonstrates metastable coexistence near critical thresholds.

## Contact

For questions or collaboration, please contact:
Jingyao Sun
jysun@igsnrr.ac.cn

# CMOEA-DRA

Source code for the manuscript *A Constrained Multiobjective Evolutionary Algorithm
with Dynamic Resource Allocation for Problems with Large Infeasible Regions*.

## Requirements
- MATLAB R2022b
- PlatEMO 4.16 — https://github.com/BIMK/PlatEMO

## Installation
Copy the `.m` files into

```
PlatEMO/Algorithms/Multi-objective optimization/CMOEA-DRA/
```

then start `platemo`.

## Run

```matlab
platemo('problem',@MW1,'algorithm',@CMOEADRA,'N',100,'maxFE',1e5,'save',1);
```

Setting used in the paper: `N = 100`, `1e5` function evaluations, 30 independent runs.

## Parameters

| # | Name | Default | Meaning |
|---|---|---|---|
| 1 | `e1` | 1 | environmental selection of FA (1 SPEA2, 2 NSGA-II, 3 modified NSGA-III) |
| 2 | `e2` | 1 | environmental selection of FEA (same coding) |
| 3 | `adapt` | 1 | enable the soft-constraint annealing in DA |
| 4 | `faHi` | 50 | FA parent quota at the start of the run |
| 5 | `faFloor` | 30 | FA parent quota at the end of the run |
| 6 | `earlyFrac` | 0.5 | fraction of the budget during which the one-way merge is active |
| 7 | `feasGate` | 25 | merge runs only while FEA holds fewer than this many feasible members |

## Notes
- `CMOEADRA` implements the algorithm called CMOEA-DRA in the paper.
- The three-archive framework is adapted from CMOEA-CD (Z. Liu et al., IEEE Trans.
  Evol. Comput. 29(6) (2025) 2771–2784). The files keep the PlatEMO copyright notice
  and are distributed for research use only.

## Citation
If you use this code, please cite the manuscript above and PlatEMO:

> Y. Tian, R. Cheng, X. Zhang, Y. Jin, PlatEMO: A MATLAB platform for evolutionary
> multi-objective optimization, IEEE Comput. Intell. Mag. 12 (4) (2017) 73–87.

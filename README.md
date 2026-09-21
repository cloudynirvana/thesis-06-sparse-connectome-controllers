# Sparse Connectome-Style Controllers as In-Silico Policy Classes

**Thesis #6** — working computational-research manuscript.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**Affiliation:** Independent computational research / Project Confluence  
**Date:** 21 September 2026

## Non-claims

This thesis is **research only**. It is not a medical device, not clinical decision support, not a protocol, not a dose, and not a cure. Sparse Kenyon-cell-style units are a cartoon of expansion recoding on a toy ODE. They are not fly neurons treating a tumour.

`fly-brain-vs-tumor` is a **game**. `malecns-immune-sight` is a **visualization**. Neither is this paper.

No DOI is registered for this document; do not invent one.

See [DISCLAIMER.md](DISCLAIMER.md).

## Problem (research-scoped)

Relative to lumped adaptive-therapy controllers, do sparse connectome-style policies change closed-loop computational behaviour on a toy cancer ordinary-differential-equation (ODE) plant in a way that is identifiable from controller architecture alone?

## Files

| Path | Role |
|---|---|
| `THESIS.md` | Full manuscript (Vancouver citations; Problem / Justification / Significance before Methods) |
| `THESIS.pdf` | Regenerated citeable PDF |
| `CITATION.cff` | Citation metadata (no document DOI) |
| `DISCLAIMER.md` | Research-only disclaimer |
| `sim/compare_controllers.py` | Toy-plant comparison (seed 20260921) |
| `sim/results.json` | Numbers reported in the manuscript |

## How to cite

Ogbonna KE. Sparse connectome-style controllers as in-silico policy classes: identifiable closed-loop differences from lumped adaptive therapy on a toy cancer ODE [Internet]. Thesis #6 working manuscript. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-06-sparse-connectome-controllers

Prefer `CITATION.cff` for machine-readable citation. When a document DOI is later minted, add it there only after it exists.

## Reproduce the toy comparison

```bash
python3 sim/compare_controllers.py
```

Requires NumPy. Output: `sim/results.json`. Trajectories are computational signatures, not patient outcomes.

## Related objects (not this paper)

| Object | Honest label |
|---|---|
| [fly-brain-vs-tumor](https://github.com/cloudynirvana/fly-brain-vs-tumor) | Browser game |
| [malecns-immune-sight](https://github.com/cloudynirvana/malecns-immune-sight) | Research visualization |
| [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub) | Catalog (this topic is NP-06) |
| [thesis-01-confluence-onco](https://github.com/cloudynirvana/thesis-01-confluence-onco) | Thesis #1 (knowledge gates) |

## Licence

Manuscript text and comparison code in this repository are provided for scholarly reuse with attribution (MIT).

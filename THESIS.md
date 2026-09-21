# Sparse Connectome-Style Controllers as In-Silico Policy Classes: Identifiable Closed-Loop Differences from Lumped Adaptive Therapy on a Toy Cancer ODE

**Document type:** Thesis #6 — working computational-research manuscript  
**Author:** Kelechi Emeka Ogbonna  
**Affiliation:** Independent computational research / Project Confluence (GitHub `cloudynirvana`)  
**Correspondence:** Kelechi Emeka Ogbonna, kelechiogbonna300@gmail.com  
**Date:** 21 September 2026  
**Status:** Architectural and methodological findings on a toy ordinary-differential-equation (ODE) plant. Not a clinical result.  
**Citation style:** numbered Vancouver. Journal DOIs were checked against Crossref on 21 September 2026. GitHub items are URL citations without a DOI.  
**DOI:** none registered for this document. Do not invent one.

This manuscript answers the control-policy identification question catalogued as NP-06 in the research theses hub. It does **not** promote the browser game `fly-brain-vs-tumor` or the visualization `malecns-immune-sight` into a paper, a protocol, or a therapy.

---

## Abstract

Adaptive-therapy controllers in mathematical oncology are typically lumped maps of a scalar burden: treat when a one-dimensional threshold is crossed, halt when it falls, and rely on competitive suppression among clones [1–9]. Expansion-recoding circuits — Kenyon-cell-style sparse codes, cerebellar granule expansions, echo-state and liquid-state reservoirs — are a different architecture class: a high-dimensional sparse feature map with a trained readout [15–32]. Those circuits are models of computation. They are not models of medical treatment.

The question of this thesis is computational, not clinical. Relative to lumped adaptive-therapy controllers, do sparse connectome-style policies change closed-loop behaviour on a toy cancer ODE in a way that is identifiable from controller architecture alone?

Three policy classes are closed around the same two-clone competitive logistic plant: (i) a Gatenby-style treat-and-halt map of total burden; (ii) an affine map of the sensitive and resistant coordinates; (iii) a frozen random-claw expansion, anterior-paired-lateral-style top-*k* sparsening, and ridge-trained readout. On a matched-burden isoline (*T* = 0.60), the lumped controller is constant (range 0; Pearson correlation of *u* with sensitive fraction 0). The sparse policy is not (range 0.91; Pearson 0.71). Opposite compositions at identical *T* produce the same lumped action (*u* = 1) and different sparse actions (*u* = 0.70 versus 0.19; Hamming distance 6 of 8 active units). Retuning the lumped hysteresis band to the sparse closed-loop burden does not recover the sparse input (root-mean-square error in *u* ≈ 0.91). Those differences are architectural signatures on a toy plant. They are not efficacy, not a dose, and not a claim that fly neurons, games, or user interfaces treat tumours.

Scientific success here is an identification result: a controller class is distinguishable from lumped adaptive therapy by closed-loop observables that cannot be absorbed into the lumped class’s parameterisation. Failure would have been a sparse policy whose input–output map is a function of *T* alone. That failure did not occur in the reported seed. The result remains a toy computational finding [33–42].

---

## Keywords

in-silico control; adaptive therapy; sparse coding; Kenyon-cell-style architecture; reservoir computing; identifiability; toy cancer ODE; connectome-inspired policy class; research-only; not a medical device

---

## Introduction

Mathematical oncology writes tumour growth, clonal competition, and treatment as dynamical systems [10–13]. Adaptive therapy, as introduced by Gatenby and colleagues, replaces continuous maximum-tolerated pulsing with a treat-for-stability heuristic: a chemosensitive majority is allowed to persist so that it competitively suppresses a less-fit resistant minority [1,4]. Subsequent work has developed multi-drug schedules, spatial agent-based tests, competitive-release models, and a list of open mathematical questions [2,3,5–9]. Optimal-control treatments of chemotherapy are older still and already treat dose as an input *u*(*t*) to an ODE [14]. That literature supplies the **lumped** controller class used here: a low-dimensional feedback map, typically of total burden, sometimes with hysteresis.

A separate literature studies sparse expansion recoding. Insect mushroom-body Kenyon cells transform a dense projection-neuron code into a sparse, decorrelated population code; feedback inhibition from the anterior paired lateral (APL) neuron is one documented sparsening mechanism [15–17,26]. Sparse coding as an efficient-representation principle is older than that circuit [18,19]. Cerebellar granule expansions [25,24], liquid-state machines [20], echo-state networks [21,23,30–32], and random-feature reservoirs [22] are the corresponding machine-learning objects: a frozen high-dimensional recurrence or expansion, a sparse or saturating nonlinearity, and a trained linear readout. Litwin-Kumar and colleagues showed that random sparse connectivity can maximise the dimension of a Kenyon-cell-like representation [27]; Babadi and Sompolinsky analysed when expansion and sparseness help a downstream classifier [28]. None of those papers is a cancer protocol.

The failure mode this thesis is written to block is a category error that is already easy to make in software: take a mushroom-body cartoon, close it on a tumour ODE, put a scoreboard on a web page, and read the result as treatment. Two public repositories adjacent to this work are exactly that kind of artefact, and they are **not this paper**. `fly-brain-vs-tumor` is a browser **game** on a toy Runge–Kutta plant [49]. `malecns-immune-sight` is a downsampled **visualization** in which a sparse Kenyon-cell code opens a toy immune ODE [50]. FlyWire and the 2024 whole-brain *Drosophila* model are connectomic and computational-neuroscience resources [43–46]. They are not tumour simulators.

The honest object is narrower. Hold the plant fixed. Vary only controller architecture. Ask whether the closed-loop input and the closed-loop state trajectory carry a signature of that architecture that a lumped adaptive-therapy parameterisation cannot absorb [33–40]. That is an identifiability question about policies, not a claim about patients.

---

## 1. Problem Statement

Relative to lumped adaptive-therapy controllers, do sparse connectome-style policies change closed-loop computational behaviour on a toy cancer ordinary-differential-equation (ODE) plant in a way that is identifiable from controller architecture alone?

That is the research problem. It is a question about in-silico controller classes. It is not a question about treating, dosing, or curing anyone [1,4,47].

Two architecture families are easy to confuse because both emit a scalar *u*(*t*) ∈ [0, 1] into the same right-hand side. A lumped adaptive-therapy controller is, in the sense used here, a hybrid map of a low-dimensional observation — typically total burden *T* = *S* + *R*, possibly with hysteresis — whose rationale in the evolutionary-therapy literature is competitive suppression, not pattern separation [1–9]. A sparse connectome-style policy is an expansion-recoding map: an observation is lifted by frozen random (or connectome-inspired) weights into a high-dimensional population, sparsened by a top-*k* or inhibitory-feedback rule, and read out by a trained linear combination [15–32]. The plant does not know which map produced *u*. The identification problem is whether an observer of (*x*(*t*), *u*(*t*)) can tell.

The problem is therefore to (i) write both families as policy classes around one toy competitive ODE; (ii) name closed-loop observables that would be identical if the sparse policy were secretly a function of *T* alone; (iii) test those observables; and (iv) refuse every interpretation that reads Kenyon cells, FlyWire neurons, games, or visualizations as medical treatment [43–46,49,50]. A sparse policy that is indistinguishable from retuned lumped hysteresis has not earned a separate architectural claim. A sparse policy that is distinguishable on a toy plant has not earned a clinic.

---

## 2. Justification of the Study

The study is justified by a documented mismatch between three literatures that already exist and a fourth object that does not.

**Lumped adaptive therapy is a specified controller class.** Gatenby, Silva, Gillies and Frieden introduced adaptive therapy as modulation to a stable tumour population rather than as maximum cell kill [1]. Zhang, Cunningham, Brown and Gatenby integrated evolutionary dynamics into a metastatic castrate-resistant prostate-cancer modelling and trial-research programme [2]. West and colleagues extended the idea to multi-drug frequency dynamics [3]. Preclinical and spatial models, competitive-release analyses, and a 2023 survey of open mathematical questions make the lumped class a real, citable design space [5–9]. That class is the baseline this thesis is obliged to beat *as an identifier*, not as a therapy.

**Sparse expansion recoding is a specified computational class.** Perez-Orive, Turner, Lin, Aso and others described sparsening and mushroom-body logic as neuroscience [15–17,26]. Olshausen and Field stated the sparse-coding principle [18,19]. Marr’s cerebellar cortex, Yamazaki and Tanaka’s cerebellar liquid-state machine, Maass’s liquid-state machine, Jaeger’s echo-state network, and the reservoir-computing unification are the control-relevant abstractions [20–25,30–32]. Litwin-Kumar et al. and Babadi and Sompolinsky give conditions under which sparse random expansion changes what a readout can separate [27,28]. Those papers justify treating “Kenyon-cell-style” as an architecture named by its computational operations (random claws, population sparseness, trained readout), not as a transplanted fly brain.

**Identifiability is the right success criterion.** Structural identifiability asks whether a unique parameter vector is consistent with noise-free input–output data [33,34,36,38]. Practical identifiability asks whether finite data actually constrain those parameters [35,40]. Observability of nonlinear biological systems is the companion question [39]. Villaverde and Banga located reverse-engineering in systems biology as a strategy with limits [37]. Kitano’s systems-biology programme is cited only as the demand that a model be a formal object, not a slogan [41,42]. Applying that literature to *controllers* rather than to plant parameters is the justification for this thesis: if two policy classes produce the same closed-loop map, architecture is not identifiable and should not be advertised.

**Existing software does not answer the question.** `fly-brain-vs-tumor` is labelled a game [49]. `malecns-immune-sight` is labelled a visualization [50]. FlyWire reconstructs an adult fly brain; Shiu et al. simulate sensorimotor processing in that brain [43–46]. None of those artefacts is a controlled comparison of policy classes against lumped adaptive therapy on a named plant with an identification test. Promoting them to a thesis would be the category error this study is justified in refusing.

The study is **not** justified as a medical device programme, a scheduling trial, a FlyWire-to-clinic translation, or a claim that sparse coding is a better cancer therapy. Those would be different objects, with different evidence.

---

## 3. Significance of the Study

The significance of this work is architectural and epistemic, not clinical.

**For mathematical oncology and control.** Adaptive-therapy papers already treat treatment as a closed-loop input [1–9,14]. They do not, as a class, ask whether a high-dimensional sparse feature map is a *different identifiable object* from a burden threshold. This thesis gives a small, reproducible test that can fail: if *u* is constant on every isoline of *T*, the sparse policy has collapsed to the lumped class. In the reported seed it did not. That is a method for keeping controller claims typed.

**For connectome-inspired machine learning.** Reservoir and mushroom-body models are justified as computers [20–32]. Closing them on a tumour ODE without an identification protocol invites a literary reading (“the fly brain treats the tumour”). The significance of stating Kenyon-cell-style operations as a policy class is that the neuroscience citations remain citations of computation. FlyWire remains a wiring diagram of a fly [43,45,46]. It does not become an oncologist.

**For this author’s public artefacts.** Theses 01–03 in this series concern knowledge gates, guideline constraints, and research objects [51–53]. They do not identify controller architecture. The game and the visualization [49,50] remain software. The significance of depositing this manuscript in its own repository is to keep that distinction inspectable: NP-06 is a paper; C-FLY and C-SIGHT are not.

**What this significance is not.** It is not clinical decision support, not a medical device, not a dosing table, not a cure, and not a claim that sparse policies control tumours better than lumped adaptive therapy. In the toy closed loop reported below, the affine composition-aware controller and the sparse controller differ from lumped treat-and-halt in *identifiable* ways; they are not ranked here by any clinical endpoint, because there is no patient and no endpoint. Translation is out of scope.

---

## Specific aims

1. Write lumped adaptive-therapy, affine composition-aware, and sparse Kenyon-cell-style maps as controller architecture classes around one toy two-clone ODE [1,15,20,21].
2. Name identification observables that would be identical if the sparse policy were secretly a function of total burden: isoline range of *u*, correlation of *u* with clone fraction at fixed *T*, immediate action at matched-burden opposite compositions, and residual error after retuning lumped hysteresis [33–40].
3. Report those observables for one frozen random expansion and one training seed, with code in `sim/compare_controllers.py`.
4. Keep FlyWire, MaleCNS-style visualizations, and the browser game outside the claim set [43–46,49,50].
5. Define success as identifiable architectural difference on a toy plant, not as disease control [47].

---

## What this thesis is

An in-silico comparison of policy classes. The plant is a competitive logistic pair (*S*, *R*) with a scalar input *u*. The sparse controller is a random-claw expansion of an eight-dimensional observation, a top-*k* sparsening, and a ridge readout trained on a composition-aware *in-silico* target. The comparison is seed `20260921`. Results are numbers from that script, not patient outcomes.

## What this thesis is not

This is not a medical device. It is not clinical decision support. It is not a protocol, a dose, a schedule, or a cure [47]. It is not FlyWire [43,45,46]. It is not the 2024 *Drosophila* computational brain model [44]. It is not `fly-brain-vs-tumor` [49]. It is not `malecns-immune-sight` [50]. It does not interpret fly neurons as treatment. It does not claim that a sparse policy is therapeutically superior to lumped adaptive therapy. It does not analyse patient data.

---

## Background

### Lumped adaptive therapy as a controller

Gatenby et al. stated the design goal as enforcing a stable tumour burden by leaving a chemosensitive population in place as a competitive suppressor of resistant cells [1]. The associated closed-loop idea is treat-and-halt on a burden threshold, not a high-dimensional feature map. Zhang et al. used a three-species evolutionary game with on/off cycles informed by tumour dynamics [2]. West et al. analysed frequency-dependent cycles and multi-drug absorbing regions [3]. Gatenby and Brown reviewed evolutionary dynamics as a framing for therapy research [4]. Enriquez-Navas et al. and Gallaher et al. tested evolutionary scheduling in preclinical and spatial computational models [5,6]. Hansen, Woods and Read treated resistance as a constraint on how a chemotherapeutic agent can be used in a model [8]. West, Ma and Newton modelled competitive release [9]. West et al. (2023) listed open questions at the mathematics–translation boundary [7]. Those papers are the lumped-class literature. They are cited as controller-design context. They are not this repository’s protocol and not a patient result.

Martin’s 1992 optimal-control scheduling paper is cited as the older control-theoretic statement that chemotherapy can be an input to an ODE [14]. Altrock, Liu and Michor, and Anderson and Quaranta, are the mathematical-oncology setting [10,11]. Nowell, and Merlo et al., are the evolutionary framing of tumour cell populations [12,13]. None of those citations identifies a parameter in the toy plant below.

### Sparse expansion as a controller architecture

In the locust and in *Drosophila*, Kenyon cells of the mushroom body receive dense olfactory projection-neuron input and emit a sparse code; oscillatory coincidence and intrinsic properties contribute to sparsening [16,17]. Lin et al. showed that Kenyon-cell–APL feedback maintains sparseness and supports discrimination of similar odours [15]. Aso et al. described mushroom-body architecture as a logic for associative learning, with Kenyon cells converging onto mushroom-body output neurons [26]. Those are neuroscience results. In this thesis they name **operations**: expansion, random-or-fixed claws, population sparseness, readout.

The same operations appear in computational theories that do not mention flies. Marr’s cerebellar cortex is an expansion onto granule cells [25]. Yamazaki and Tanaka treated the cerebellum as a liquid-state machine [24]. Maass, Natschläger and Markram defined liquid-state machines as real-time computers without stable states [20]. Jaeger and Haas, and the echo-state literature, defined reservoirs with trained readouts [21,30–32]. Lukoševičius and Jaeger, and Verstraeten et al., unified those methods as reservoir computing [22,23]. Pathak et al. showed that reservoirs can forecast spatiotemporal chaos as a computational task [32]. Sussillo and Abbott trained chaotic recurrent nets to produce coherent outputs [29]. Olshausen and Field located sparseness as a coding principle [18,19]. Litwin-Kumar et al. and Babadi and Sompolinsky supplied the representation-dimension and classification analysis for sparse expansions [27,28].

The architectural claim allowed here is only this: those operations define a policy class. The architectural claim forbidden here is that an insect neuron, a FlyWire edge, or a game sprite is a therapeutic actuator.

### Connectomes and why they are not this plant

Dorkenwald et al. published a neuronal wiring diagram of an adult fly brain and the FlyWire community resource [43,45]. Schlegel et al. annotated cell types across that connectome [46]. Shiu et al. built a leaky-integrate-and-fire model of *Drosophila* sensorimotor processing from connectivity and neurotransmitter identity [44]. Those papers are cited so that a later reader cannot pretend this thesis used them as a tumour engine. It did not. No FlyWire weight is an ODE coefficient here. No MaleCNS reconstruction is simulated here. The sparse controller below has 96 units and three random claws per unit. That is a cartoon of expansion recoding, comparable in spirit to the stylized controller in the game repository [49], and it is **not** a downsampled FlyWire.

### Identifiability of policies

Bellman and Åström defined structural identifiability [33]. Ljung and Glad treated global identifiability for arbitrary parametrizations [34]. Raue et al. separated structural from practical identifiability via profile likelihood [35]. Chis, Banga and Balsa-Canto compared methods on systems-biology models [36]. Villaverde and colleagues reviewed reverse engineering, structural identifiability, and observability [37–39]. Wieland et al. restated the distinction for contemporary practice [40]. In this manuscript the “parameters” of interest are not plant rate constants. They are the policy class. The lumped class is parametrized by hysteresis thresholds (*T*<sub>on</sub>, *T*<sub>off</sub>). The affine class is parametrized by a hyperplane in (*S*, *R*). The sparse class is parametrized by a frozen expansion and a readout. Two classes are indistinguishable, for the purpose of this thesis, if every closed-loop observable of the sparse policy can be reproduced by some parameter of the lumped class. Distinguishing them is the result. Fitting the plant to data is not attempted [35,40].

---

## Methods

### Honesty gates

| Gate | Requirement |
|---|---|
| 0 | Research-only scope. Not a medical device. Not clinical decision support [47]. |
| 1 | Plant is a toy ODE. Symbols *S*, *R*, *T* are unitless computational burdens, not tumour volumes in a person. |
| 2 | Controller classes are named by their maps, not by clinical intent. |
| 3 | Fly neurons, games, and visualizations are not treatment [49,50]. |
| 4 | Identification observables are functions of (*x*, *u*), not of survival, RECIST, or dose [33–40]. |
| 5 | Numbers come from `sim/compare_controllers.py` with seed 20260921, or are refused. |
| 6 | No document DOI is minted. Journal DOIs are Crossref-checked. GitHub items have no `doi:` field. |

### Plant

The plant is a two-clone competitive logistic system with a single bounded input *u* ∈ [0, 1]:

\[
\frac{dS}{dt} = r_S S \left(1 - \frac{S + \alpha_{SR} R}{K}\right) - \delta_S u\, S
\]

\[
\frac{dR}{dt} = r_R R \left(1 - \frac{R + \alpha_{RS} S}{K}\right) - \delta_R u\, R
\]

with \(T = S + R\). Parameters used in the reported seed are \(r_S = 0.28\), \(r_R = 0.16\), \(\alpha_{SR} = 1.0\), \(\alpha_{RS} = 1.6\), \(\delta_S = 0.55\), \(\delta_R = 0.06\), \(K = 1\). Time is an arbitrary computational unit (“days” only as an integrator label). Integration is classical fourth-order Runge–Kutta with step 0.05 on \([0, 200]\). The cost-of-resistance pattern (\(r_R < r_S\), \(\delta_R \ll \delta_S\), \(\alpha_{RS} > \alpha_{SR}\)) is a qualitative cartoon of the competitive-suppression story in the adaptive-therapy literature [1,8,9]. It is not identified from data. It is not a named human tumour.

### Controller class A — lumped adaptive therapy

A hybrid treat-and-halt map of total burden [1]:

- if \(T \ge T_{\mathrm{on}} = 0.50\), set \(u = 1\);
- if \(T \le T_{\mathrm{off}} = 0.25\), set \(u = 0\);
- otherwise hold the previous \(u\).

This is the lumped class. It cannot, by construction, distinguish two states with the same \(T\) and different \((S, R)\) except through subsequent dynamics. After a reset, \(u = 0\), so open-loop evaluation on an isoline inside the hysteresis band reports \(u = 0\) everywhere on that isoline.

### Controller class B — affine composition-aware map

\[
u = \mathrm{clip}(w_S S + w_R R + b,\, 0,\, 1)
\]

with \((w_S, w_R, b) = (2.4, -1.8, -0.55)\). This class sees composition and does **not** expand. It is the control for the confound “any policy that reads *S* and *R* will differ from lumped AT.” If the sparse class were only that confound, it would be redundant with class B.

### Controller class C — sparse Kenyon-cell-style policy

Observation (eight projection-neuron-style channels):

\[
\mathrm{PN} = \bigl(S,\, R,\, T,\, 1-T,\, S/T,\, R/T,\, u_{\mathrm{prev}},\, 1\bigr).
\]

A population of \(N = 96\) units each draws \(n_{\mathrm{claw}} = 3\) channels without replacement, with independent Gaussian weights (mean 0, variance 1) and a small Gaussian bias. Drive is the claw inner product. The top \(k = 8\) units (APL-style population sparseness [15]) keep a rectified drive; the rest are zero; the code is \(\ell_2\)-normalised. The readout is

\[
u = \mathrm{clip}(h^\top a + a_0,\, 0,\, 1),
\]

with \(a\) fit by ridge regression (\(\lambda = 10^{-2}\)) on 800 random states in \((0.02, 0.85)^2\). The training *target* is an in-silico composition-aware function, not the lumped AT map: treat when \(T \ge 0.50\) and the sensitive fraction is at least 0.35; stay off when \(T \le 0.22\) or the sensitive fraction is below 0.20. That target is a computational object used to give the readout something composition-dependent to implement. It is not a clinical rule and is not copied from a trial [1,2].

Weights of the expansion are frozen after sampling (reservoir-style [21–23]). Only the readout is trained. Seed 20260921.

### Identification protocol

The protocol is designed so that a negative result is possible.

1. **Isoline test.** For fixed \(T \in \{0.40, 0.60\}\), vary the sensitive fraction on \((0.05, 0.95)\) and record \(u\) after controller reset. A policy that is a function of \(T\) alone has range 0 and Pearson correlation 0 against sensitive fraction.
2. **Matched-burden opposite compositions.** Immediate \(u\) at \((S, R) = (0.48, 0.12)\) versus \((0.12, 0.48)\), both \(T = 0.60\), and at the more extreme pair \((0.50, 0.10)\) versus \((0.10, 0.50)\). Hamming distance of the active Kenyon-cell-style set is reported for the latter pair.
3. **Decision-region summary.** On a 41 × 41 grid of \((S, R) \in [0.02, 0.90]^2\), the fraction of states with \(u > 0.5\) and the number of 4-connected ON components.
4. **Closed loop.** Three initial conditions: matched-*T* sensitive-rich, matched-*T* resistant-rich, and a mid mix. Report duty cycle of \(u\), switch count of the Boolean \(u > 0.5\), mean \(T\), and terminal resistant fraction. These are dynamical signatures, not efficacy.
5. **Retuning transfer.** Fit lumped \((T_{\mathrm{on}}, T_{\mathrm{off}})\) by grid search to a one-dimensional reduction of the sparse closed-loop burden on the sensitive-rich trajectory, then simulate the retuned lumped controller on both matched-*T* initial conditions. Large residual RMSE in \(u\) and in \(T\) means the lumped parameterisation cannot absorb the sparse architecture.

### What was not done

No patient data. No FlyWire or MaleCNS weights. No pharmacokinetics identified from plasma. No toxicity constraint calibrated to a host. No comparison against a fitted 15-dimensional CONFLUENCE plant [51]. No claim that the sparse readout is optimal. No hyperparameter sweep beyond the reported seed. No medical device validation.

---

## Results

These are computational identification findings on a toy plant. They are not experimental oncology results and not patient outcomes [47]. All numbers are from `sim/compare_controllers.py`, seed 20260921.

### The lumped map is constant on burden isolines; the sparse map is not

At \(T = 0.60\), lumped adaptive therapy emits \(u = 1\) for every composition (range 0, standard deviation 0, Pearson correlation of \(u\) with sensitive fraction 0). The affine map has range 0.76 and Pearson 0.78. The sparse Kenyon-cell-style map has range 0.91 (minimum 0.06, maximum 0.97, standard deviation 0.30) and Pearson 0.71. At \(T = 0.40\), which lies inside the lumped hysteresis band after reset, the lumped map is identically 0; the sparse map still varies (range 0.64).

On the \(T = 0.60\) isoline, the sparse policy disagrees with the lumped ON/OFF threshold at 33% of sampled compositions (mean absolute difference in \(u\) of 0.42). The affine map disagrees at 88% of those points. The sparse class is therefore not a copy of lumped treat-and-halt, and it is not a copy of the affine hyperplane either: decision-region ON fractions are 0.86 (lumped), 0.18 (affine), and 0.73 (sparse). Each ON set was a single 4-connected component in this grid; component count is **not** claimed as a discriminator here.

### Opposite compositions at identical burden produce identical lumped actions and different sparse actions

Immediate actions at \(T = 0.60\):

| Initial condition | \((S, R)\) | Lumped AT | Affine | Sparse KC-style |
|---|---|---|---|---|
| Sensitive-rich | (0.48, 0.12) | 1.00 | 0.39 | 0.70 |
| Resistant-rich | (0.12, 0.48) | 1.00 | 0.00 | 0.19 |
| Mid mix (\(T = 0.50\)) | (0.30, 0.20) | 1.00 | 0.00 | 0.57 |

At the opposite-composition pair with the same \(T = 0.60\), \((0.50, 0.10)\) versus \((0.10, 0.50)\), lumped \(\Delta u = 0\), affine \(\Delta u = 0.47\), sparse \(\Delta u = 0.53\). The sparse active sets differ in 6 of 8 units (Hamming distance 6). That is a pattern-separation signature of the expansion, in the limited sense of this cartoon [15,27,28]. It is not a claim about odour discrimination in a fly.

### Closed-loop signatures differ; they are not ranked as therapies

Over 200 time units from the sensitive-rich matched-*T* initial condition, lumped treat-and-halt has duty cycle 0.97, two ON/OFF switches, mean \(T = 0.57\), and terminal resistant fraction 1.00 (sensitive clone numerically extinguished). The affine map has duty cycle 0.18, no Boolean switches, mean \(T = 0.64\), and terminal resistant fraction 0.31. The sparse map has duty cycle 0.10, three switches, mean \(T = 0.83\), and terminal resistant fraction 0.88. From the resistant-rich matched-*T* initial condition the lumped duty cycle is 1.00 with \(\Delta u = 0\) at \(t = 0\) relative to the sensitive-rich case; the sparse duty cycle is 0.051 with immediate \(u = 0.19\) rather than 0.70. From the mid-mix initial condition the sparse controller switches seven times; the affine controller switches zero times.

These trajectories show that architecture changes the closed-loop signal. They do **not** show that any class “controls cancer better.” Mean burden is higher under the sparse policy in this seed. A reader who converts that sentence into a dosing recommendation has left the paper.

### Lumped hysteresis cannot be retuned to absorb the sparse input

Grid search of \((T_{\mathrm{on}}, T_{\mathrm{off}})\) against a one-dimensional reduction of the sparse sensitive-rich burden path returned the extreme pair \((0.80, 0.05)\). Closed-loop RMSE against the sparse policy was 0.28 in \(T\) and 0.91 in \(u\) on that initial condition, and 0.29 in \(T\) and 0.94 in \(u\) on the resistant-rich transfer condition. Duty cycles remained near 1 for the retuned lumped controller versus 0.10 and 0.05 for the sparse controller. The lumped class, as parametrized here, does not contain the sparse class.

That is the identification result. It is local to this plant, this seed, and this hysteresis parametrization. It is not a general theorem [33,34,38].

---

## Discussion

The sparse policy differs from lumped adaptive therapy on every observable that would have been identical if *u* were a function of *T* alone. That is what the problem statement asked. The affine controller also differs from lumped adaptive therapy, as it must, because it is a function of \((S, R)\). The sparse controller is not redundant with the affine controller: isoline range, ON fraction, closed-loop switch count, and Hamming change under opposite compositions are different. Architecture, in this toy, is not a synonym for “uses two clones.”

The result should be read with the same coldness as a failed identifiability test in the other direction. If the ridge readout had been trained to imitate lumped treat-and-halt, the isoline range would have been near zero and the architectural claim would have been false. The training target was deliberately composition-aware so that the expansion had something non-lumped to implement. That is a feature of the methods, and it is also a limitation: the thesis shows that a sparse expansion *can* implement an identifiable non-lumped policy, not that every sparse expansion does. A randomly wired Kenyon-cell cartoon with a readout trained on \(T\) alone would be a lumped controller in expensive clothing.

Relation to the cited literatures is correspondingly narrow. Gatenby’s treat-for-stability idea explains why a lumped burden threshold is a serious baseline rather than a straw man [1,4]. It does not license this sparse map as adaptive therapy. Reservoir computing explains why a frozen expansion plus a trained readout is a legitimate computer [20–23]. It does not license this plant as a physical reservoir of a tumour. Lin et al. and Litwin-Kumar et al. explain why sparseness and claw count are not arbitrary decorations [15,27]. They do not explain cancer.

`fly-brain-vs-tumor` and `malecns-immune-sight` remain a game and a visualization [49,50]. A later citation that treats those repositories as clinical evidence, or that treats this thesis as their validation, is a mis-citation. FlyWire and Shiu et al. remain fly resources [43–46]. This manuscript does not use their weights, cell types, or spike trains.

Theses 01–03 in this series refuse skip-level promotion of knowledge into parameters [51–53]. This thesis refuses a different skip: promotion of a neural motif into a treatment, and promotion of a game into a paper. The conversion ladder is the same honesty rule with a different first object.

---

## Limitations

- One plant, one seed, one claw sample, one readout training draw [47].
- The lumped class is a two-threshold hysteresis, not the full design space of multi-drug adaptive therapy [3,7].
- The sparse class is 96 random-claw units, not FlyWire, not MaleCNS, and not a leaky-integrate-and-fire mushroom body [43,44].
- No measurement noise model on *S* and *R*; an observer who sees only *T* cannot implement classes B or C as written [39].
- No structural-identifiability proof for the readout; the argument is empirical on named observables [33,35,38,40].
- Closed-loop “performance” is not an endpoint. Higher or lower mean *T* in this seed is not ranked.
- Adjacent software (game, visualization) is not re-analysed here and is not a result [49,50].
- No wet-lab measurement, no patient, no trial.

---

## Future work

1. Repeat the isoline and retuning tests across expansion seeds and claw counts, including a readout trained to imitate lumped AT as an explicit negative control [27].
2. Give all controllers the same observation sigma-algebra (for example *T* only) and ask whether expansion of a scalar still produces identifiable closed-loop differences, or whether composition sensing is doing all the work.
3. Replace hysteresis with a broader lumped class (dose modulation, vacation-oriented rules [6]) before claiming architecture against “adaptive therapy” as a whole [7].
4. Keep FlyWire weights out of the plant unless a later paper states a *computational* hypothesis that a named connectome statistic changes a named identification observable — and still without a treatment claim [43–46].
5. Optional later: mint a document DOI after a preprint deposit; do not invent one here.

Translation to care is out of scope [47].

---

## Conclusions

On a toy two-clone ODE, a sparse Kenyon-cell-style policy is not a lumped adaptive-therapy controller. Isoline range, composition-conditioned action at matched burden, active-set Hamming distance, closed-loop duty cycle, and residual error after hysteresis retuning all distinguish the architecture classes in seed 20260921. An affine map of \((S, R)\) is a third class, not a synonym for either.

That is an in-silico identification result. It is not evidence that fly neurons, games, visualizations, or connectomes treat tumours. `fly-brain-vs-tumor` and `malecns-immune-sight` are not this paper. Adaptive-therapy and reservoir-computing literatures remain what they were: evolutionary scheduling as a lumped control idea [1–9], and sparse expansion as a computer [15–32]. Holding those objects apart is the point of the thesis.

---

## References

1. Gatenby RA, Silva AS, Gillies RJ, Frieden BR. Adaptive therapy. Cancer Res. 2009;69(11):4894-4903. doi:10.1158/0008-5472.CAN-08-3658.
2. Zhang J, Cunningham JJ, Brown JS, Gatenby RA. Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. Nat Commun. 2017;8:1816. doi:10.1038/s41467-017-01968-5.
3. West J, You L, Zhang J, Gatenby RA, Brown JS, Newton PK, et al. Towards multidrug adaptive therapy. Cancer Res. 2020;80(7):1578-1589. doi:10.1158/0008-5472.CAN-19-2669.
4. Gatenby RA, Brown JS. Integrating evolutionary dynamics into cancer therapy. Nat Rev Clin Oncol. 2020;17(11):675-686. doi:10.1038/s41571-020-0411-1.
5. Enriquez-Navas PM, Kam Y, Das T, Hassan S, Silva A, Foroutan P, et al. Exploiting evolutionary principles to prolong tumor control in preclinical models of breast cancer. Sci Transl Med. 2016;8(327):327ra24. doi:10.1126/scitranslmed.aad7842.
6. Gallaher JA, Enriquez-Navas PM, Luddy KA, Gatenby RA, Anderson ARA. Spatial heterogeneity and evolutionary dynamics modulate time to recurrence in continuous and adaptive cancer therapies. Cancer Res. 2018;78(8):2127-2139. doi:10.1158/0008-5472.CAN-17-2649.
7. West J, Adler F, Gallaher J, Strobl M, Brady-Nicholls R, Brown J, et al. A survey of open questions in adaptive therapy: bridging mathematics and clinical translation. Elife. 2023;12:e84263. doi:10.7554/eLife.84263.
8. Hansen E, Woods RJ, Read AF. How to use a chemotherapeutic agent when resistance to it threatens the patient. PLoS Biol. 2017;15(2):e2001110. doi:10.1371/journal.pbio.2001110.
9. West J, Ma Y, Newton PK. Capitalizing on competition: an evolutionary model of competitive release in metastatic castration resistant prostate cancer treatment. J Theor Biol. 2018;455:249-260. doi:10.1016/j.jtbi.2018.07.028.
10. Altrock PM, Liu LL, Michor F. The mathematics of cancer: integrating quantitative models. Nat Rev Cancer. 2015;15(12):730-745. doi:10.1038/nrc4029.
11. Anderson ARA, Quaranta V. Integrative mathematical oncology. Nat Rev Cancer. 2008;8(3):227-234. doi:10.1038/nrc2329.
12. Merlo LMF, Pepper JW, Reid BJ, Maley CC. Cancer as an evolutionary and ecological process. Nat Rev Cancer. 2006;6(12):924-935. doi:10.1038/nrc2013.
13. Nowell PC. The clonal evolution of tumor cell populations. Science. 1976;194(4260):23-28. doi:10.1126/science.959840.
14. Martin R. Optimal control drug scheduling of cancer chemotherapy. Automatica. 1992;28(6):1113-1123. doi:10.1016/0005-1098(92)90054-J.
15. Lin AC, Bygrave AM, de Calignon A, Lee T, Miesenböck G. Sparse, decorrelated odor coding in the mushroom body enhances learned odor discrimination. Nat Neurosci. 2014;17(4):559-568. doi:10.1038/nn.3660.
16. Perez-Orive J, Mazor O, Turner GC, Cassenaer S, Wilson RI, Laurent G. Oscillations and sparsening of odor representations in the mushroom body. Science. 2002;297(5580):359-365. doi:10.1126/science.1070502.
17. Turner GC, Bazhenov M, Laurent G. Olfactory representations by Drosophila mushroom body neurons. J Neurophysiol. 2008;99(2):734-746. doi:10.1152/jn.01283.2007.
18. Olshausen BA, Field DJ. Emergence of simple-cell receptive field properties by learning a sparse code for natural images. Nature. 1996;381(6583):607-609. doi:10.1038/381607a0.
19. Olshausen BA, Field DJ. Sparse coding of sensory inputs. Curr Opin Neurobiol. 2004;14(4):481-487. doi:10.1016/j.conb.2004.07.007.
20. Maass W, Natschläger T, Markram H. Real-time computing without stable states: a new framework for neural computation based on perturbations. Neural Comput. 2002;14(11):2531-2560. doi:10.1162/089976602760407955.
21. Jaeger H, Haas H. Harnessing nonlinearity: predicting chaotic systems and saving energy in wireless communication. Science. 2004;304(5667):78-80. doi:10.1126/science.1091277.
22. Lukoševičius M, Jaeger H. Reservoir computing approaches to recurrent neural network training. Comput Sci Rev. 2009;3(3):127-149. doi:10.1016/j.cosrev.2009.03.005.
23. Verstraeten D, Schrauwen B, D'Haene M, Stroobandt D. An experimental unification of reservoir computing methods. Neural Netw. 2007;20(3):391-403. doi:10.1016/j.neunet.2007.04.003.
24. Yamazaki T, Tanaka S. The cerebellum as a liquid state machine. Neural Netw. 2007;20(3):290-297. doi:10.1016/j.neunet.2007.04.004.
25. Marr D. A theory of cerebellar cortex. J Physiol. 1969;202(2):437-470. doi:10.1113/jphysiol.1969.sp008820.
26. Aso Y, Hattori D, Yu Y, Johnston RM, Iyer NA, Ngo TT, et al. The neuronal architecture of the mushroom body provides a logic for associative learning. Elife. 2014;3:e04577. doi:10.7554/eLife.04577.
27. Litwin-Kumar A, Harris KD, Axel R, Sompolinsky H, Abbott LF. Optimal degrees of synaptic connectivity. Neuron. 2017;93(5):1153-1164.e7. doi:10.1016/j.neuron.2017.01.030.
28. Babadi B, Sompolinsky H. Sparseness and expansion in sensory representations. Neuron. 2014;83(5):1213-1226. doi:10.1016/j.neuron.2014.07.035.
29. Sussillo D, Abbott LF. Generating coherent patterns of activity from chaotic neural networks. Neuron. 2009;63(4):544-557. doi:10.1016/j.neuron.2009.07.018.
30. Jaeger H. Echo state network. Scholarpedia. 2007;2(9):2330. doi:10.4249/scholarpedia.2330.
31. Yildiz IB, Jaeger H, Kiebel SJ. Re-visiting the echo state property. Neural Netw. 2012;35:1-9. doi:10.1016/j.neunet.2012.07.005.
32. Pathak J, Hunt B, Girvan M, Lu Z, Ott E. Model-free prediction of large spatiotemporally chaotic systems from data: a reservoir computing approach. Phys Rev Lett. 2018;120(2):024102. doi:10.1103/PhysRevLett.120.024102.
33. Bellman R, Åström KJ. On structural identifiability. Math Biosci. 1970;7(3-4):329-339. doi:10.1016/0025-5564(70)90132-X.
34. Ljung L, Glad T. On global identifiability for arbitrary model parametrizations. Automatica. 1994;30(2):265-276. doi:10.1016/0005-1098(94)90029-9.
35. Raue A, Kreutz C, Maiwald T, Bachmann J, Schilling M, Klingmüller U, et al. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. Bioinformatics. 2009;25(15):1923-1929. doi:10.1093/bioinformatics/btp358.
36. Chis OT, Banga JR, Balsa-Canto E. Structural identifiability of systems biology models: a critical comparison of methods. PLoS One. 2011;6(11):e27755. doi:10.1371/journal.pone.0027755.
37. Villaverde AF, Banga JR. Reverse engineering and identification in systems biology: strategies, perspectives and challenges. J R Soc Interface. 2014;11(91):20130505. doi:10.1098/rsif.2013.0505.
38. Villaverde AF, Barreiro A, Papachristodoulou A. Structural identifiability of dynamic systems biology models. PLoS Comput Biol. 2016;12(10):e1005153. doi:10.1371/journal.pcbi.1005153.
39. Villaverde AF. Observability and structural identifiability of nonlinear biological systems. Complexity. 2019;2019:8497093. doi:10.1155/2019/8497093.
40. Wieland FG, Hauber AL, Rosenblatt M, Tönsing C, Timmer J. On structural and practical identifiability. Curr Opin Syst Biol. 2021;25:60-69. doi:10.1016/j.coisb.2021.03.005.
41. Kitano H. Systems biology: a brief overview. Science. 2002;295(5560):1662-1664. doi:10.1126/science.1069492.
42. Kitano H. Computational systems biology. Nature. 2002;420(6912):206-210. doi:10.1038/nature01254.
43. Dorkenwald S, Matsliah A, Sterling AR, Schlegel P, Yu SC, McKellar CE, et al. Neuronal wiring diagram of an adult brain. Nature. 2024;634(8032):124-138. doi:10.1038/s41586-024-07558-y.
44. Shiu PK, Sterne GR, Spiller N, Franconville R, Sandoval A, Zhou J, et al. A Drosophila computational brain model reveals sensorimotor processing. Nature. 2024;634(8032):210-219. doi:10.1038/s41586-024-07763-9.
45. Dorkenwald S, McKellar CE, Macrina T, Kemnitz N, Lee K, Lu R, et al. FlyWire: online community for whole-brain connectomics. Nat Methods. 2022;19(1):119-128. doi:10.1038/s41592-021-01330-0.
46. Schlegel P, Yin Y, Bates AS, Dorkenwald S, Eichler K, Brooks P, et al. Whole-brain annotation and multi-connectome cell typing of Drosophila. Nature. 2024;634(8032):139-152. doi:10.1038/s41586-024-07686-5.
47. Ogbonna KE. DISCLAIMER.md [Internet]. thesis-06-sparse-connectome-controllers / GitHub; 2026 Sep 21 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-06-sparse-connectome-controllers/blob/main/DISCLAIMER.md
48. Ozturk MC, Xu D, Príncipe JC. Analysis and design of echo state networks. Neural Comput. 2007;19(1):111-138. doi:10.1162/neco.2007.19.1.111.
49. Ogbonna KE. fly-brain-vs-tumor [Internet]. GitHub; 2026 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/fly-brain-vs-tumor
50. Ogbonna KE. malecns-immune-sight [Internet]. GitHub; 2026 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/malecns-immune-sight
51. Ogbonna KE. CONFLUENCE × OnCo: an evidence-gated dynamical framework for integrating oncology knowledge graphs with adaptive cancer-state models [Internet]. Thesis #1 working manuscript. 2026 Sep 20 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-01-confluence-onco
52. Ogbonna KE. Complexity science and NSTG-guided in-silico pathology dynamics for biologics pathway exploration [Internet]. Thesis #2 computational research manuscript. 2026 Sep 20 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-02-complexity-nstg
53. Ogbonna KE. Disease profiles for complex pathologies: a gated method for systemic personalized-medicine research objects [Internet]. Thesis #3 working method manuscript. 2026 Sep 20 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-03-disease-profile
54. Ogbonna KE. Research theses hub [Internet]. GitHub; 2026 Sep 21 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/research-theses-hub
55. Ogbonna KE. sim/compare_controllers.py [Internet]. thesis-06-sparse-connectome-controllers / GitHub; 2026 Sep 21 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-06-sparse-connectome-controllers/blob/main/sim/compare_controllers.py

Journal DOIs appear only for Crossref-verified records. GitHub items are complete Internet citations without a `doi:` field. No DOI is invented.

---

## Disclaimer

**Research technical report.** This repository is a computational research manuscript. It is not a medical device, not a clinical decision-support system, not a diagnostic or therapeutic product, and not a protocol [47]. This document makes no cure claim, no dosing recommendation, and no claim of patient benefit. Simulated trajectories are not patient outcomes. Sparse Kenyon-cell-style units in the accompanying script are a cartoon of expansion recoding. They are not fly neurons acting on a tumour. `fly-brain-vs-tumor` is a game [49]. `malecns-immune-sight` is a visualization [50]. FlyWire and related whole-brain models are neuroscience resources [43–46]. None of those objects is this paper, and none is medical treatment.

Author: Kelechi Emeka Ogbonna — kelechiogbonna300@gmail.com

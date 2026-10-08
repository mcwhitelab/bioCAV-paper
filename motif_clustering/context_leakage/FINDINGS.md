# Does whole-protein context leak into domain-window CAVs?

Run 2026-09-06. Library: **/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all**
(27,472 motifs, L25). NOT the /groups AHMAD_CAV_OUTPUTS_DIR (7,692 motifs, missing
9 of the 11 preATP-grasp families).

## 1. The pipeline does not trim before the model  (confirmed)

`esm2_tcav/src/embed.py::process_dataset` appends `(accession, sample['sequence'])`
-- the **full protein** -- and runs it through ESM-C. `pool_window_embeddings` then
slices `hidden_states[start+1 : end+1]`. Trimming is **pooling-only**; self-attention
saw the whole chain. So context leakage has a live mechanism.

Background set: `neg_data/neg_10000.jsonl`, embedded once to
`_shared_negatives/L25_all.npy` (10,000 x 1152, float64) and **symlinked** as
`L25_neg.npy` into every motif dir. All 27,459 CAVs share one identical negative
cloud. Those negatives are pooled over the **whole protein**
(`batch_train_cavs.py::prepare_negative_jsonl_full_sequence` sets
`window_span = [0, sequence_length]`), while positives are pooled over a
~100-200 aa domain window.

## 2. The dendrogram observation is real

`cut_k00314` cluster 289 has 54 of 27,459 leaves: 29 ATP-grasp (CL0179),
17 NADP_Rossmann (CL0063), 8 unassigned. It contains **13/13 preATP-grasp
families and 28/43 ATP-grasp families**. (Median k=314 cluster size is 3,
max 9,737 -- so 289 is not the grab-bag cluster.)

At `cut_k02769` it splits into a near-pure preATP-grasp cluster (267, n=20:
every preATP-grasp domain in Pfam, spanning CL0063 + CL0179 + unassigned)
and an ATP-grasp cluster (2182, n=18).

## 3. But raw CAV cosine contradicts the simple leakage story

|                              | mean cos |
|------------------------------|---------:|
| within preATP-grasp (n=11)   |  +0.3473 |
| within ATP-grasp (n=43)      |  +0.1275 |
| within other CL0063 (n=307)  |  +0.0544 |
| **preATP <-> ATP-grasp**     |**-0.0353**|
| preATP <-> other CL0063      |  +0.0514 |
| preATP <-> 800 random        |  +0.0339 |
| ATP-grasp <-> 800 random     |  +0.0317 |

preATP-grasp CAVs are **not** aligned with ATP-grasp CAVs -- their cosine is at
or below the random baseline. What is real is that preATP-grasp forms a very
tight group of its own (+0.347, far above its own clan's +0.054).

Cognate same-protein pairs are likewise *dis*similar, while ATP-grasp families in
*different* proteins are similar (all figures recomputed on the 27k xdisk library;
the two libraries agree wherever both carry the family):

    cognate (same protein)              control (different proteins)
    Dala_Dala_lig_N | _C     -0.556     ATP-grasp   | ATP-grasp_3   +0.782
    Rimk_N          | RimK   -0.496     Dala_Dala_lig_C | RimK      +0.737
    GARS_N          | GARS_A -0.179     GARS_A      | GSH-S_ATP     +0.353
    Biotin_carb_N | CPSase_L_D2 -0.139
    Ins134_P3_kin_N | core   +0.064
    CPSase_L_D1     | D2     +0.039
    Synapsin      | Synapsin_C +0.181

Mean cognate -0.155 vs mean control +0.624. If host-protein context dominated the
window embedding, cognate pairs -- which share a host protein -- should be the
*similar* ones. They are the dissimilar ones.

So the k=314 merge is not driven by CAV direction similarity. The dendrogram uses
diffusion (walktrap) distance on a mutual-kNN graph, where groups can be close via
shared neighbours without direct cosine similarity.

## 4. Three live explanations

1. **Context leakage** (original hypothesis). Mechanism confirmed in (1), but the
   direct cosine evidence in (3) does not support it in its simple form.
2. **Genuine remote homology.** Structurally the ATP-grasp fold is an N-terminal
   (= "pre-ATP-grasp") subdomain plus the ATP-grasp core; Pfam's CL0063/CL0179
   split is a Pfam convention. ESM may be recognising a real superfamily.
   *Needs checking against ECOD/SCOP before it can be ruled in or out.*
3. **Diffusion-graph artifact** of leaf ordering / coarse merges.

## 5. Next experiment (supersedes the cognate-pair cosine design)

cos(CAV_pre, CAV_grasp) is the wrong readout -- it is ~0 to begin with, so it has
no room to collapse. Use instead a direct, CAV-free measure of leakage:

**Host-architecture decodability.** For one Pfam family, split its positives by
host-protein architecture (e.g. host also contains ATP-grasp vs not). Train a
classifier to predict host architecture *from the domain-window embedding alone*.
- full-context regime: AUROC >> 0.5 means context leaks, quantitatively.
- excised regime (domain substring alone): AUROC should fall to 0.5.

This measures leakage magnitude directly, needs one family and ~500 sequences,
and does not depend on the preATP/ATP-grasp story being true.

## Caveat on `run_excision.py` as currently written

It gives each family its own negative windows (length-matched per family). That
breaks the library's constant-background design, which is what makes CAV cosines
comparable across motifs -- it injects variance into the very quantity being
measured. If that script is revived, use a single fixed negative window set shared
by all families within each regime.

---

# 6. PF07065 / PF25596 case study (2026-09-06)

## 6a. A natural experiment that settles the leakage question

Carbamoyl-phosphate synthase large chain carries **PF25596** (CPSase_L_D1,
preATP-grasp, clan CL0063) immediately N-terminal to **PF02786** (CPSase_L_D2,
ATP-grasp, clan CL0179). **PF00289** (Biotin_carb_N) is the *same* preATP-grasp
domain type but in a *different* host family (biotin / pyruvate / urea
carboxylase, not CPS). That gives a clean 2x2:

|                                    | cos    |
|------------------------------------|-------:|
| CPSase_L_D1 x CPSase_L_D2  (**same protein**, different domain type) | **+0.039** |
| CPSase_L_D1 x Biotin_carb_N (**different protein**, same domain type) | **+0.519** |
| Biotin_carb_N x CPSase_L_D2 (different protein, different domain type) | -0.139 |

Domain identity beats host-protein identity by >13x. If whole-protein context
dominated the window embedding, the first row would be the large one. It is not.

Top-10 library neighbours confirm it -- each family's neighbourhood is its own
structural type, drawn from many different host proteins:
- PF25596: Biotin_carb_N, GARS_N, RS_preATP-grasp-like, PylC-like_N, ATPgrasp_N,
  DUF1246, Dala_Dala_lig_N -- all preATP-grasp.
- PF00289: GARS_N +0.708, RS_preATP-grasp-like +0.631, PylC-like_N +0.620,
  Dala_Dala_lig_N +0.558 ... all preATP-grasp.
- PF02786: GARS_A +0.593, ATP-grasp +0.543, ATP-grasp_3 +0.482, ATP-grasp_4,
  RimK, Dala_Dala_lig_C -- all CL0179.

**Conclusion: these CAVs are domain detectors, not architecture detectors.**
This closes the context-leakage hypothesis for the preATP/ATP-grasp block.

## 6b. PF07065 (D123) is a different kind of anomaly -- an untrimmed CAV

D123's annotated domain covers a **median 83% of its host protein**
(domain 294 aa, protein 348 aa). There is essentially no context to leak: its
"domain-window" CAV is effectively a whole-protein CAV.

Its neighbour list shows the cost. The top four are real (TTL +0.455,
Ins134_P3_kin +0.454, GSH_synth_ATP +0.433, Synapsin_C +0.386 -- ATP-grasp), but
positions 5-10 are unrelated eukaryotic proteins of similar length with no clan:
DCA16 +0.353, DUF789 +0.335, COPR5 +0.328, IFTAP +0.327, Serglycin +0.327,
CCDC32 +0.312. That is generic "eukaryotic protein of ~350 aa" signal, and it is
the failure mode to expect wherever domain coverage approaches 100%.

Note Pfam itself flags D123 as ATP-grasp-derived ("ATP-grasp enzymes belonging to
the R2K clade ... part of the core ATP-grasp module"), so its ATP-grasp
neighbours are correct; the junk tail is the artifact.

## 6c. PF25596 has a positives-quality problem

**63.6% of its 1,000 positives are UniProt fragments**, median host length 320 aa
against a real CPS large chain of ~1,070. (PF00289: 3.5% fragments, median host
633 aa.) So PF25596's CAV is trained largely on truncated hosts -- a caveat on
that CAV independent of anything above, and a reason it is a **singleton** at
cut_k02769 (cluster 1805, n=1) despite cos +0.519 to PF00289.

Worth screening `%fragment` and `domain coverage of host` across all 27k
families: they are cheap per-CAV quality flags, and both are computable from the
existing `pos_1000.jsonl` files with no GPU.

---

# 7. What IS driving the block  (2026-09-06)

> **SUPERSEDED by sections 8a and 9.** The "N-terminal subdomain" framing below is
> wrong -- 8a shows nine CL0179 families share that exact geometry without leaning
> pre. The block is a fold signal, and the merge with ATP-grasp is a seam effect
> (section 9), not a positional one. Section 7a (not a diffusion artifact) still
> stands. Read 7b as descriptive only, and treat 7d as retired.

## 7a. Not a diffusion artifact

The same block appears under **cosine+average** (leaf ranks 1263-1286) as under
diffusion (4336-4357), with a near-identical membership list and the same
terminating neighbour (PF18086 PPIP5K2_N). Method is not the cause.

(The earlier "median gap" statistic -- 0.2555 for cosine_average vs 0.0007 for
diffusion -- was misleading. ATP-grasp occupies 25,703 of 27,459 leaves, so its
*median* leaf position is meaningless. Discard that comparison.)

## 7b. The block is a pure pre-ATP-grasp / N-terminal-subdomain family
##     that crosses Pfam clan boundaries

Median domain start position within the host protein:

    IN THE BLOCK                                   CONTROLS (ATP-grasp core)
    PF02951 GSH-S_N        CL0063    0.6%  N-term  PF02222 ATP-grasp     29.7% central
    PF00289 Biotin_carb_N  CL0063    0.0%  N-term  PF07478 Dala_Dala_C   38.5% central
    PF02844 GARS_N         CL0063    0.0%  N-term  PF08443 RimK          32.3% central
    PF21360 PylC-like_N    CL0063    1.6%  N-term
    PF27100 DcsG_N         -         1.7%  N-term
    PF03738 GSP_synth      CL0179    2.9%  N-term   <-- ATP-grasp clan!
    PF03917 GSH_synth_ATP  CL0063    3.3%  N-term
    PF08886 GshA           -         3.7%  N-term
    PF18086 PPIP5K2_N      CL0063    3.9%  N-term
    PF18130 ATPgrasp_N     -        10.4%  N-term
    PF14403 CP_ATPgrasp_2  CL0179   15.0%  N-term   <-- ATP-grasp clan!
    (exceptions: PF18419 ATP-grasp_6 56.6% C-term; PF03199 GSH_synthase 43.5% central)

12 of 15 measurable members start in the N-terminal fifth of their host. The block
spans **three different Pfam clan labels** -- CL0063, CL0179 and unassigned -- and
the CL0179 members in it are ones Pfam itself describes as pre-ATP-grasp
(PF03738 = "Glutathionylspermidine synthase preATP-grasp"; CP_ATPgrasp_1/2 =
circularly permuted ATP-grasp).

**So the CAVs recover a coherent structural family that Pfam's clan boundaries
split.** That is the finding -- ESM being consistent where Pfam is not.

## 7c. The "clusters together with ATP-grasp" part is a coarse-cut effect

preATP core vs ATP-grasp core mean cosine = **-0.202** (each is internally +0.48
and +0.43). The best bridging family is PF03738 GSP_synth at min(to_pre,to_grasp)
= +0.093, about 1 SD above the global background of +0.036 +/- 0.056. Nothing
pulls them together.

They land in one 54-member cluster only at `cut_k00314` (289), where two compact
but mutually distant groups become each other's least-bad merge partner. At
`cut_k02769` they are cleanly separate (267 vs 2182). The circle-pack adjacency is
what catches the eye; it is not similarity.

## 7d. RETIRED -- positional leakage (rejected in section 8a)

*Kept for the record. Do not test this; section 8a falsified it.*

Host-protein *identity* does not drive these CAVs (section 6a: same-protein
different-domain +0.039 vs different-protein same-domain +0.519). But **position
in the chain** might: this block is overwhelmingly N-terminal domains, and
N-terminal position is exactly the kind of information full-context embedding
supplies and excision would destroy.

This was proposed as testable and then falsified in section 8a before being run:
nine CL0179 families (RLAN, RNA_lig_T4_1, DNA_ligase_aden, RNA_ligase, PPDK_N,
Citrate_synth_N, ATP-grasp_2, CPSase_L_D2, GSH_synth_ATP) share the N-terminal,
short geometry and do not lean pre. Domain start-fraction does not explain the
grouping. **Do not run this test.**

---

# 8. Do "combined" pre+ATP-grasp families cluster with pure pre-ATP-grasp?  (2026-09-06)

Tested directly. Families split by domain geometry (median over their 1,000
positives), scored as **lean = mean cos to preATP core - mean cos to ATP-grasp
core**. Data: `combined_vs_split.csv`.

| group | n | mean lean | median lean |
|---|--:|--:|--:|
| A. COMBINED (domain >= 60% of host) | 8 | **-0.222** | -0.099 |
| B. N-terminal & short (pure preATP-like geometry) | 23 | +0.386 | +0.592 |
| C. central / C-terminal | 23 | -0.351 | -0.455 |

**The claim is not supported as stated.** Combined-domain families lean toward
ATP-grasp, not pre-ATP-grasp: ATPgrasp_YheCD -0.640, ATPgrasp_ST -0.544,
ATPgrasp_TupA -0.425, D123 -0.167.

**But the underlying observation is right, for a different reason.** Four of the
eight -- CP_ATPgrasp_2 +0.056, GSP_synth +0.002, GshA -0.025, DUF3182 -0.030 --
have lean ~ 0 *and* weak similarity to **both** cores (to_pre +0.05..+0.10,
to_grasp +0.03..+0.09, against a global background of +0.036 +/- 0.056). They are
ambivalent, not pre-leaning. That is exactly why they sit at the *boundary* of the
preATP block in leaf order (diffusion ranks 4333-4337) rather than inside either
core.

Mechanism: mean-pooling a window that spans both subdomains averages the two
signals toward each other and cancels them. A combined-domain CAV is a **diluted**
vector, not a vector reassigned to one subdomain. Near-zero lean != preATP-leaning.

## 8a. This retires the positional hypothesis of section 7d

Group B contains nine CL0179 families with the *same* N-terminal, short geometry
as the pre-ATP-grasp families -- RLAN, RNA_lig_T4_1, DNA_ligase_aden, RNA_ligase,
PPDK_N, Citrate_synth_N, ATP-grasp_2, CPSase_L_D2, GSH_synth_ATP -- and they do
**not** lean pre (lean 0 to -0.52). N-terminal position alone does not produce the
pre-ATP-grasp grouping. **Section 7d's positional-leakage hypothesis is wrong and
should not be tested further.** The preATP block is a genuine fold signal.

## 8b. A concrete falsifiable prediction

Exactly **two** families have a CAV lean that disagrees with their Pfam clan:

- **PF18419 (ATP-grasp_6)** -- filed CL0179, but to_pre **+0.384**, to_grasp
  **-0.190** (lean +0.574), and it sits inside the preATP cluster 267.
- **PF03917 (GSH_synth_ATP)** -- filed CL0063 (Rossmann/preATP), but it is the
  ATP-binding half of eukaryotic glutathione synthetase; lean +0.031, ambivalent
  (see section 9a).

Not misassigned, despite earlier drafts saying so: PF03199 (lean +0.370) and
PF25596 (lean +0.341) are both already in CL0063, which is where their lean puts
them. They are consistent with Pfam, not corrections to it.

A structural check (ECOD/SCOP, or Foldseek against a known pre-ATP-grasp
subdomain) on **PF18419** alone would settle whether this is a genuine
"the CAVs corrected Pfam" result. Optional follow-up, not a blocker.

---

# 9. What joins the two cores: misannotated + merged families  (2026-09-06)

Rebuilt the pipeline's own graph (`cav_cluster.knn_graph`, k=15, floor=0.10,
mutual) over all 27,459 CAVs: 88,946 undirected edges, mean degree 6.5,
186 components. Both cores land in component 109.

**Shortest kNN path, preATP core -> ATP-grasp core = 4 hops, and every
intermediate node is a misannotated or merged family:**

    [PRE core]    PF02078  Synapsin        CL0063
      seam        PF03199  GSH_synthase    CL0063   misannotation
      seam        PF03917  GSH_synth_ATP   CL0063   misannotation
      seam        PF07065  D123            CL0179   merged (83% of host)
    [GRASP core]  PF05770  Ins134_P3_kin   CL0179

## 9a. The eukaryotic glutathione synthetase case

PF03199 ("Eukaryotic glutathione synthase") and PF03917 ("Eukaryotic glutathione
synthase, ATP binding domain") are the two halves of one protein -- the
pre-ATP-grasp half and the ATP-grasp half -- and **Pfam files both in CL0063**
(NADP_Rossmann). Compare the prokaryotic enzyme, which Pfam splits correctly:
PF02951 GSH-S_N -> CL0063, PF02955 GSH-S_ATP -> CL0179.

The CAVs detect the inconsistency: PF03199 lean **+0.370** (pre), PF03917 lean
**+0.031** (ambivalent -- it is the ATP-grasp half wearing a Rossmann label).
Host proteins confirm the pairing: both are ~99% "Glutathione synthetase",
PF03917 spanning 3-44% of the host and PF03199 44-65%.

## 9b. The seam families (leaf ranks 4332-4339, contiguous)

    PF03199 GSH_synthase   CL0063  cov 21%  lean +0.370  misannotated clan
    PF03917 GSH_synth_ATP  CL0063  cov 40%  lean +0.031  misannotated clan
    PF14403 CP_ATPgrasp_2  CL0179  cov 76%  lean +0.064  merged + circ. permuted
    PF28486 AnkG_N         CL0179  cov 48%  lean +0.166  (geometry from
                                                          _prepared_jsonl, n=900)
    PF08886 GshA           -       cov 94%  lean -0.024  merged
    PF03738 GSP_synth      CL0179  cov 96%  lean +0.011  merged; Pfam names it
                                                          "...preATP-grasp"
    PF04174 CP_ATPgrasp_1  CL0179  cov  ?   lean +0.152  circularly permuted --
                              NO positives on disk (absent from tcav_data and from
                              _prepared_jsonl); classification rests on the Pfam
                              name alone, not on measured geometry
    PF25596 CPSase_L_D1    CL0063  cov 25%  lean +0.341  correct; k2769 singleton
    PF18419 ATP-grasp_6    CL0179  cov  7%  lean +0.574  filed CL0179, CAV says pre

Two mechanisms, exactly as hypothesised by the user:
1. **Merged domains** -- one Pfam model spans both subdomains, so mean-pooling
   averages the two signals and yields a diluted, ambivalent CAV (section 8).
   Circular permutation (CP_ATPgrasp_1/2) is a variant of the same problem.
2. **Clan misannotation** -- families whose CAV lean disagrees with their Pfam
   clan (PF03917, PF18419, and arguably PF03738 whose own Pfam name says
   preATP-grasp while its clan says ATP-grasp).

## 9c. Consequence for the paper

The k=314 merge of cluster 267 with 2182 is not a similarity statement (the cores
are anti-correlated at -0.202). It is produced by a thin chain of ambiguous
families. That makes the observation a **feature, not a bug**: the CAV graph
localises exactly those Pfam entries whose domain boundaries or clan assignments
are questionable, and the seam is short enough to enumerate by hand.

## 9d. The seam is load-bearing — quantified

Removing families and re-measuring the shortest kNN path between the two cores:

    full graph                          4 hops
    drop PF03199+PF03917 (euk GS pair)  7 hops   <-- load-bearing
    drop PF07065 (D123, merged)         4 hops   <-- parallel route exists
    drop all 3 path nodes               7 hops
    drop all 10 seam families           7 hops

The cores are never fully disconnected, so the seam is not the *only* link. But
scale matters. In component 109 (26,881 nodes) the median random-pair distance is
**9 hops**; only **1.39%** of random pairs are within 4 hops, while **20.7%** are
within 7.

So the 4-hop seam route is a **1st-percentile** connection -- genuinely
exceptional. Remove it and the two cores fall back to 7 hops, which is ordinary
background proximity. And the residual 7-hop route is a generic ramble through
unrelated nucleotide-binding families with no ATP-grasp content:

    RS_preATP-grasp-like -> DXP_reductoisom -> UDPG_MGDP_dh_N -> Mannitol_dh
      -> Lar_N -> PmbA_TldD_2nd -> S-AdoMet_synt_M -> GARS_A

**Conclusion.** The only *meaningful* connection between the pre-ATP-grasp and
ATP-grasp clusters runs through Pfam families that are either clan-misannotated
(PF03917, PF18419, PF03738) or whose domain model merges both subdomains
(PF07065, PF08886, PF14403, PF04174, PF03738). Strip those and the two groups are
no closer than any random pair of nucleotide-binding families.

---

# 10. Final assessment  (2026-09-06, after structural confirmation of PF18419)

**Structural check confirmed PF18419 (ATP-grasp_6) is pre-ATP-grasp, not
ATP-grasp** -- agreeing with its CAV (lean +0.574) against its Pfam clan (CL0179).

## 10a. Cluster 267 is 100% pre-ATP-grasp

PF18419 was the *only* CL0179-labelled member of cluster 267 (n=20; the other 19
are CL0063 or unassigned). With the structural correction the cluster has **no
ATP-grasp members at all**. The premise of the original observation -- that
pre-ATP-grasp and ATP-grasp families cluster *together* -- does not survive: at
k=2769 the pre-ATP-grasp cluster is pure.

## 10b. What actually draws the two groups together: nothing

Cores are **anti-correlated** (-0.202) while internally cohesive (+0.483, +0.427).
There is no attraction. Three separate things created the appearance:

1. **Leaf-order adjacency is a rendering property.** In a 27,459-leaf dendrogram,
   any two clusters that merge anywhere must be drawn side by side. This is what
   is visible in the circle-pack, and it is not a similarity statement.
2. **The k=314 merge** is produced by a thin seam of ambiguous Pfam entries --
   a 4-hop kNN path (1st percentile; median random pair 9 hops) through
   PF03199 -> PF03917 -> PF07065. Strip the seam and the cores sit 7 hops apart,
   which is ordinary background (20.7% of random pairs). Two failure modes:
   clan misannotation (PF03917) and merged/circularly-permuted domain models
   (PF07065, PF08886, PF14403, PF04174, PF03738), whose mean-pooled CAVs are
   *diluted* and therefore weakly near both cores.
3. **One apparent contaminant was Pfam's error, not contamination** (10a).

Ruled out along the way, each with a measurement:
- **Context leakage** (the original hypothesis): same-protein/different-domain
  cos +0.039 vs different-protein/same-domain +0.519 (section 6a).
- **Positional / N-terminal leakage**: nine CL0179 families share the geometry and
  do not lean pre (section 8a).
- **Method artifact**: the block reproduces under cosine+average (section 7a).
- **Remote homology between the clans**: the cores are anti-correlated, and the
  only short path between them runs through annotation errors (section 9).

## 10c. The claim to make in the paper

Invert it. This is not "CAV clustering is contaminated by protein context." It is:

> **CAV clustering recovers the pre-ATP-grasp fold as a coherent group spanning
> three Pfam clan labels, and the few places where it disagrees with Pfam are
> Pfam's errors** -- PF18419 (confirmed structurally), PF03917 (the ATP-grasp half
> of eukaryotic glutathione synthetase, filed in the Rossmann clan while the
> prokaryotic enzyme is split correctly), and merged domain models whose CAVs
> come out diluted.

The seam is short enough to enumerate by hand, which makes this a usable audit
signal, not a caveat: **short kNN paths between otherwise anti-correlated CAV
clusters localise questionable Pfam domain boundaries and clan assignments.**

# FamFlow — Protein Family Exploration

FamFlow takes a researcher from a **Seed** to a set of **Community** descriptions and, optionally, profile HMMs. It is a chain of independent **Stages** that communicate only through files and one **Sequence Table**, not a fixed pipeline. Biological judgement (what to include, which cutoff, which markers count) stays with the researcher; FamFlow runs the computation and reports what it finds.

This is a single bounded context (`famflow.py`). ROTIFER, HMMER/pyhmmer, mmseqs, and NCBI are upstream systems (see *Neighbouring systems*).

## Language

### Workflow

**Stage**:
One subcommand (`seed`, `search`, … `scan`). It reads one defined artifact, writes one defined artifact plus a **Stage Report**, and knows nothing about later stages.
_Avoid_: step, phase, module (the code comments say "estagio")

**Terminal Stage**:
The principle that every **Stage** is a valid place to stop, because each one leaves behind a readable **Stage Report**.
_Avoid_: checkpoint (that means something else here)

**Stage Report**:
The human-readable `REPORT.md` a **Stage** writes into its output directory.
_Avoid_: log, summary

**Run Record**:
The history of the exact parameters of every run of a **Stage** (`.<stage>_params.json`). The **Workflow Report** reads it.
_Avoid_: config, settings

**Checkpoint**:
Per-**Stage** record of which items (usually **Units**) are done and which failed, with tracebacks, so that re-running a command resumes the work (`.<stage>_checkpoint.json`, the `Progress` class).
_Avoid_: progress, state

**Workflow Report**:
The single Markdown document (`report` stage) combining Methods (from **Run Records**) and Results (from **Stage Reports**).
_Avoid_: relatorio (in code), final report

### Inputs and search

**Unit**:
The named thing being pursued (a domain, a protein, or a whole family). Every per-file artifact is keyed by it (`<unit>.fasta`, `<unit>.tsv`, …).
_Avoid_: domain, query, target, dataset (the README says "domain", but the code deliberately stays neutral)

**Seed**:
The starting material for a **Unit**: a single sequence, a set of sequences, a ready-made HMM, or one sequence cut into **Domains** by a **Chopping**.
_Avoid_: query, input

**Chopping**:
A spec that splits one seed sequence into discontinuous **Domains**, for example `1-169_323-361,171-318` (commas separate domains, underscores separate segments).
_Avoid_: split, segmentation

**Domain**:
One **Unit** made by a **Chopping** (`<base>_dominio_<i>`). The README uses the word more loosely.

**Search**:
An HMM search of a **Unit**'s profile against the **Reference Database**. It reports **Hits** at a deliberately permissive **Reporting E-value**.

**Reference Database**:
The protein database that **Search** runs against and that sequences are fetched from (default: the NCBI nr path in `NR`).
_Avoid_: db, bank ("banco")

**Hit**:
One row of **Search** output: a target sequence aligned to a **Unit**'s profile, with score, e-value, and envelope coordinates.
_Avoid_: match, result

**Profile Coverage**:
The fraction of the *HMM profile* covered by a **Hit** (`qcov`). This is the quantity that `--min-coverage` filters on. It is not **Target Coverage**.
_Avoid_: coverage (without a qualifier), cov

**Cutoff Sensitivity**:
The grid showing how many **Hits** survive each combination of **Profile Coverage** and e-value. `inspect` produces it before any filter is applied.

### Collection

**Collected Sequence**:
A **Hit** that passed the `collect` filters and was fetched from the **Reference Database**. Collected Sequences become rows of the **Sequence Table**.
_Avoid_: hit (after collection), member

**Collection Mode**:
Either **Envelope** or **Full**. This choice determines whether the network compares a fold or whole architectures.

**Envelope**:
Only the aligned region of a target, plus an `--expand` margin. The network then compares the same fold.
_Avoid_: domain region, crop

**Full**:
The entire target protein. Use it when **Communities** might differ in domain architecture.

**Target Length / Target Coverage**:
The length of the whole source protein, and the fraction of it that a **Collected Sequence** covers. A low Target Coverage means the sequence is an isolated domain inside a larger multidomain protein.
_Avoid_: talilen (that is the aligned length, not the protein length)

**Anchor**:
The **Seed** sequence, guaranteed to be present and identifiable in a **Unit**'s collected set, so that its **Community** can be highlighted later. It is optional, and no **Community** is privileged without one.
_Avoid_: query, reference (but see *Flagged ambiguities*)

**Anchor Resolution**:
How the **Anchor** entered the set: *exact self-hit* (renamed from an identical **Hit**), *assumed* (the longest **Hit** at or above `--anchor-identity`), or *inserted* (added as a new sequence).

### Redundancy and network

**Cluster**:
A group of near-identical **Collected Sequences** from mmseqs redundancy reduction at a given identity and coverage. It is a bookkeeping device, not a biological group.
_Avoid_: community, group, family

**Representative**:
The one **Collected Sequence** that stands for its **Cluster**. Only Representatives go on to `matrix` and `ssn`.
_Avoid_: centroid, seed

**Cluster Size**:
How many **Collected Sequences** a **Representative** stands for. It is the only valid abundance measure; counting Representatives does not estimate abundance.

**Cluster Sweep**:
A diagnostic grid of identity × coverage showing how many **Clusters** each combination produces. It is run before choosing thresholds and writes no FASTA.

**Similarity Matrix**:
The all-vs-all pairwise hits between the **Representatives** of one **Unit**. Its e-value is a *ceiling* on possible edges, not the final cutoff.
_Avoid_: distance matrix, blast table

**SSN (Sequence Similarity Network)**:
A graph of one **Unit**'s **Representatives**, with an edge wherever a pair's bitscore is at or above the **Edge Cutoff**.
_Avoid_: network, graph (without a qualifier)

**Edge Cutoff**:
The bitscore threshold that creates SSN edges. By default it is chosen automatically at the peak of the **Closeness Scan**.
_Avoid_: threshold, score cutoff (that belongs to `inspect`)

**Closeness Scan**:
The curve of network closeness across candidate **Edge Cutoffs**. It can be *truncated* (still rising where the scan ends) or *flat* (no peak).

**Single-Community Fallback**:
What happens when the **Closeness Scan** has no reliable peak (N is too small): the whole **Unit** becomes one **Community**, and no edges are fabricated.

**Community**:
A partition of an **SSN** found by community detection (Louvain/Leiden). It is the primary result object: every Community gets a **Community Profile**, and none is treated as noise.
_Avoid_: cluster, family, group, clan

**Orphan**:
A **Representative** with no edge above the **Edge Cutoff**. Orphans are the most divergent candidates, not garbage.
_Avoid_: singleton (that is a **Cluster** of size 1)

**Family**:
A *hypothesis* that a **Community** is separable by a profile HMM. **Separability QC** tests it. Nothing in the code is a Family until then.

### Annotation and genomic context

**Annotation**:
Best-effort taxonomy (organism → genus/lineage) and domain **Architecture** added to the **Sequence Table**. An Annotation failure never fails the workflow.

**Architecture**:
The ordered domain composition of a protein (`arch`, `pfam`, `aravind` columns from ROTIFER).

**Neighborhood**:
The genes within `--window` positions on either side of a **Collected Sequence**'s gene, from NCBI. One Neighborhood is one `block_id`.
_Avoid_: context (without a qualifier, it collides with "bounded context"), operon, locus

**Scope**:
Which **Collected Sequences** get a **Neighborhood** lookup: `all`, `anchor` (the **Anchor**'s **Community**), `community:X`, or `top:N`.

**Marker Vocabulary**:
A researcher-supplied JSON file of named regex **Markers** (with a class: `mandatory`, `accessory`, …) and **Vetoes**. It encodes the biology that the code deliberately leaves out.
_Avoid_: profile config, rules, markers file

**Marker**:
A named **Neighborhood** gene (for example `gspD`) whose presence supports membership in a biological system.

**Veto**:
A named **Neighborhood** gene whose presence contradicts membership (for example `pilT` for T2SS).

**Context Evaluation**:
The Marker and Veto *counts* for one **Neighborhood**, without a verdict. Only **Selection** turns these counts into accept or reject.

### Description, selection, models

**Community Profile**:
One row per **Community** comparing size (**Representatives** and summed **Cluster Size**), lengths, **Target Coverage**, mean SSN degree, top taxa, top **Architectures**, and **Marker** statistics (`communities.tsv`).
_Avoid_: summary, perfil (in code)

**Discovery Mode**:
The `profile` behaviour when no **Marker Vocabulary** is given: it ranks **Enriched Terms** for each **Community** against the background, which finds the vocabulary instead of requiring one.

**Enriched Term**:
A **Neighborhood** annotation term that is over-represented in one **Community** relative to all Communities (add-1 log2 odds). It is exploratory, not a statistical test.

**Selection**:
The optional filter that turns **Context Evaluations** into ACEITO (accepted) or REJEITADO (rejected) **Decisions** and writes one **Curated Set** per **Community**. It restores the behaviour of the old query-centric pipeline.
_Avoid_: filter, curation

**Decision**:
The accept or reject verdict for one sequence, with a reason (`decisoes.tsv`).

**Curated Set**:
The accepted sequences of one **Community**, used as input to `build` (`<community>_curado.fasta`).

**Model**:
A profile HMM built from an aligned **Curated Set**. It carries its **Community**'s name.
_Avoid_: profile (which is overloaded with the `profile` stage)

**Model Library**:
All **Models** concatenated and hmmpressed into one file (default `familias.hmm`).

**Self-Scan**:
Scanning every **Model** against the pooled **Representatives** to find each sequence's best-scoring Model.

**Separability QC**:
The fraction of sequences whose best **Model** in the **Self-Scan** is their own **Community**'s Model. A low value means the **Communities** are not **Families**.

### Central artifact

**Sequence Table**:
The single table with one row per **Collected Sequence** (`00_table/sequences.tsv`), keyed by `id`. Stages only *add or overwrite columns*, merging on `id`. It is the deliverable, and the column is the contract: any **Stage** can be replaced by an external tool that fills the same column.
_Avoid_: central table ("tabela central"), metadata, dataframe

## Relationships

- A **Seed** defines exactly one **Unit**, except that a **Chopping** makes one **Unit** (a **Domain**) per chopped region.
- A **Unit** has one **Search**, which produces many **Hits**.
- `collect` turns **Hits** into **Collected Sequences** (at most one per target, the best score). This is the stage that creates rows in the **Sequence Table**.
- A **Unit** has at most one **Anchor**.
- Every **Collected Sequence** belongs to exactly one **Cluster**, and each **Cluster** has exactly one **Representative**.
- A **Unit** has one **Similarity Matrix** and one **SSN**, built over its **Representatives** only.
- An **SSN** partitions its **Representatives** into one or more **Communities**. **Orphans** stay outside the graph unless `--keep-orphans` is set.
- A **Community** has one **Community Profile**, at most one **Curated Set**, and at most one **Model**.
- A **Neighborhood** belongs to one **Collected Sequence**. When several match, the Neighborhood with the most **Markers** wins.
- **Separability QC** compares **Model** names to **Community** names in the **Sequence Table**.

### Column ownership in the Sequence Table

| Stage | Columns it writes |
|---|---|
| `collect` | `id`, `unit`, `length`, `target_length`, `target_coverage`, `search_evalue`, `search_score`, `search_qcov`, `is_anchor` (from FASTA: `organism`, `description`) |
| `cluster` | `cluster_rep`, `cluster_size`, `is_cluster_rep` |
| `ssn` | `community`, `ssn_degree`, `is_orphan` |
| `annotate` | taxonomy (`taxid`, `lineage`, `phylum`…`genus`) without overwriting existing values; `pfam`, `aravind`, `arch` |
| `dashboard` | reads `is_reference`, which nothing writes (see below) |

## Neighbouring systems

- **ROTIFER** (`rotifer.devel.alpha.*`, `rotifer.devel.beta.sequence`, `rotifer.db.ncbi`, `rotifer.genome.data`): the *Conformist* upstream. FamFlow consumes its DataFrame shapes directly (for example `pyhmmer_to_df` columns and the `c{cov}i{id}` cluster column name). The only translation is `COLUMN_ALIASES` / `normalize_search_columns` (`cov` → `coverage`). SSN building, the **Closeness Scan**, and the dashboard live in ROTIFER (`mvroliveira`), not here.
- **HMMER / pyhmmer, mmseqs, diamond/blast, famsa**: tools used through ROTIFER or `subprocess`.
- **NCBI**: the upstream source of taxonomy and **Neighborhoods**.

## Example dialogue

> **Dev:** "After `cluster`, is a **Community** just a bigger **Cluster**?"
> **Domain expert:** "No. A **Cluster** only removes redundancy: near-identical sequences collapse into one **Representative**. **Communities** come later, from the **SSN** over those Representatives. Different question, different threshold."

> **Dev:** "So to say how big a **Community** is, I count its rows?"
> **Domain expert:** "Count its rows and you get **Representatives**. Sum the **Cluster Size** to get sequences. A Community of 5 Representatives can stand for 5,000 sequences."

> **Dev:** "The **Anchor**'s Community is the answer, and the rest is noise?"
> **Domain expert:** "That was the old pipeline. Here every **Community** gets a **Community Profile**. The Anchor only tells you which one to highlight. And **Orphans** are the most interesting part, not the least."

> **Dev:** "When `profile` reports **Markers**, has it accepted those sequences?"
> **Domain expert:** "No. A **Context Evaluation** is just counts. Only **Selection** makes a **Decision**, and only if you run it. If the **Models** built from the **Curated Sets** then fail **Separability QC**, those Communities weren't **Families**. Revisit the **Edge Cutoff**."

## Flagged ambiguities

- **"coverage" has four meanings**: **Profile Coverage** (`qcov`, `--min-coverage` in `inspect`/`collect`), **Target Coverage** (`target_coverage`), mmseqs alignment coverage for a **Cluster** (`cluster --coverage`), and the raw input column `coverage`/`cov`. `inspect` also overwrites `df["coverage"]` with `qcov`. Resolution: always qualify the word, and rename `cluster --coverage` to something like `cluster_coverage`.
- **"Unit" vs "domain" vs "family"**: the README says *domain*; the code says *unidade* (unit); `seed --chopping` produces `_dominio_` units. Resolution: **Unit** is the general term, **Domain** is only a chopped Unit, and **Family** is only a tested hypothesis about a **Community**.
- **"Cluster" vs "Community"**: the `cluster` stage reports "quantas comunidades em cada agrupamento" ("how many communities in each grouping") and the sweep heatmap is titled with clusters. Resolution: these are distinct concepts (see the dialogue), and `cluster` output should never say "community".
- **Anchor vs Reference**: `collect` writes `is_anchor`, but `dashboard` looks for `is_reference`. No stage writes `is_reference`, even though a comment claims `seed` records it. Resolution: these mean the same thing; use **Anchor** and fix `dashboard`.
- **Community names are not unique across Units**: each **Unit**'s **SSN** labels its communities independently (`Louvain_0`, …; the **Single-Community Fallback** always uses `Louvain_0`). `profile`, `select`, `context --scope community:X`, and **Separability QC** all group by the bare `community` column across the whole **Sequence Table**, so communities from different Units with the same label get merged. Resolution for a redesign: **Community** identity is (**Unit**, label).
- **Orphans in the Sequence Table**: the `ssn` Stage Report says Orphans are marked `is_orphan` in the table. That only happens with `--keep-orphans`, because without it Orphans are never nodes and get no row in the merge. Resolution: decide whether `is_orphan` is always written.
- **`hits_dir` means two things**: in `cluster` it is `03_hits` (**Collected Sequences**); in `matrix`, `ssn`, `select`, and `scan` it defaults to `04_clustered` (**Representatives**). Resolution: name the parameter after the concept it carries.
- **E-value has five roles**: **Search** reporting (`--evalue`), **Search** inclusion (`--inc-evalue`), the `collect` filter (`--max-evalue`), the **Similarity Matrix** edge ceiling, and the `scan` threshold. Each should be named for its role.
- **Marker class is ignored by Selection**: the **Marker Vocabulary** declares classes (`mandatory`, `accessory`), but **Selection** counts only the total `n_markers`, so "mandatory" enforces nothing. Resolution pending: either make **Selection** class-aware or drop the class from the format.
- **Model ↔ Community naming**: `select` writes the **Curated Set** file under a sanitized **Community** name, and **Separability QC** compares that sanitized name back to the raw community label. The two match only if sanitizing changed nothing.
- **"profile"** means both the HMM profile (as in **Profile Coverage**) and the `profile` stage that produces **Community Profiles**. Resolution: use **Model** for the HMM and **Community Profile** for the stage's output.
- **"context"**: in this codebase it means a genomic **Neighborhood** (`context` stage, `--no-context`), not a DDD bounded context. Prefer **Neighborhood** in new code.

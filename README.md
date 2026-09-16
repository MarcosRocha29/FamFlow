# FamFlow

**FamFlow** is a solution for optimizing, rather than automating, the serial construction of HMM models by integrating ROTIFER functions into a common workflow for processing multiple protein domains.

The project integrates functions from the **ROTIFER** ecosystem to apply sequence-processing and alignment steps consistently across multiple protein domains. The aim is to reduce repetitive manual operations while keeping the biological decisions involved in model construction under the user's control.

## Overview

The construction of HMM models for protein domains commonly involves several repeated steps, particularly when multiple domains are being investigated.

FamFlow organizes these steps into a common workflow:

```text
Protein/domain sequences
          │
          ▼
   Sequence processing
          │
          ▼
      Alignment
          │
          ▼
 Alignment inspection
          │
          ▼
   HMM model building
```

The same workflow can be applied to multiple domains:

```text
             ┌── Domain 1 ──► processing ──► alignment ──► HMM
             │
Input data ──┼── Domain 2 ──► processing ──► alignment ──► HMM
             │
             ├── Domain 3 ──► processing ──► alignment ──► HMM
             │
             └── Domain N ──► processing ──► alignment ──► HMM
```

FamFlow does not determine whether a sequence or domain should be included in a model. Instead, it provides a framework for carrying out the computational steps associated with that process.

## Motivation

When building several HMMs, the same sequence-processing and alignment operations may need to be repeated for every domain.

FamFlow is intended to make these operations easier to reproduce and manage by providing a single workflow for multiple datasets.

The main objectives are:

* Reduce repetitive computational work
* Apply the same processing steps across domains
* Integrate existing ROTIFER functionality
* Keep intermediate files and results organized
* Facilitate reproducible HMM construction workflows

## ROTIFER integration

FamFlow uses functionality from **ROTIFER** for sequence manipulation and alignment-related operations.

Rather than implementing these operations independently, FamFlow acts as a workflow layer that coordinates them across multiple domains.

```text
                    FamFlow
                       │
          ┌────────────┴────────────┐
          │                         │
      Domain A                   Domain B
          │                         │
          ▼                         ▼
      ROTIFER                    ROTIFER
          │                         │
          ▼                         ▼
      Alignment                  Alignment
          │                         │
          └────────────┬────────────┘
                       ▼
                  HMM building
```

## Workflow

A typical analysis may include:

1. Prepare domain-specific sequence datasets.
2. Load and process the sequences.
3. Generate or manipulate multiple sequence alignments.
4. Inspect and refine the alignments.
5. Prepare alignments for HMM construction.
6. Build HMM models using the resulting alignments.
7. Compare and evaluate the resulting models.

The exact sequence of operations depends on the analysis being performed.

## Multiple sequence alignments

Multiple sequence alignment is an important intermediate step in the workflow.

FamFlow can organize alignment-related operations for multiple domains so that equivalent processing steps can be performed without manually repeating the same procedure for each dataset.

The resulting alignments can then be inspected before proceeding to HMM construction.

## HMM construction

FamFlow is intended to **optimize the workflow surrounding HMM construction**, rather than replace the researcher’s decisions about model definition.

For example, decisions concerning:

* sequence inclusion;
* domain boundaries;
* removal of problematic sequences;
* alignment quality;
* selection of representative sequences; and
* interpretation of the resulting model

remain part of the analysis.

The computational workflow can then be applied consistently once these decisions have been made.

## Sequence Similarity Networks

Sequence Similarity Networks (SSNs) can also be used as a complementary analysis when investigating protein families.

FamFlow can organize information associated with sequence-similarity analyses, allowing relationships between sequences and families to be examined alongside sequence-level information.

A general analysis can therefore combine:

```text
Sequence data
     │
     ├──────────────► Sequence similarity
     │                       │
     │                       ▼
     │                      SSN
     │                       │
     │                       ▼
     │                 Family structure
     │
     └──────────────► Multiple alignment
                             │
                             ▼
                         HMM model
```

This can be useful when investigating divergent protein families or selecting groups of sequences for domain-specific HMM construction.

## Network analysis

For SSN-based analyses, FamFlow can provide information useful for exploring network structure, including cluster relationships and network centrality measures.

These analyses are complementary to the alignment/HMM workflow and can help characterize sequence relationships within a protein family.

## Input

Depending on the workflow, FamFlow can work with:

* Protein sequences
* Domain-specific sequence datasets
* Multiple sequence alignments
* Sequence-similarity data
* Similarity matrices
* Sequence Similarity Networks
* Sequence metadata

## Output

Depending on the analysis, outputs may include:

```text
Multiple sequence alignments
Processed sequence datasets
HMM-ready alignments
HMM models
Similarity matrices
SSN/network data
Cluster information
Network statistics
```

The exact files produced depend on the command and workflow being executed.

## Installation

Clone the repository:

```bash
git clone git@github.com:MarcosRocha29/FamFlow.git
cd FamFlow
```

Alternatively:

```bash
git clone https://github.com/MarcosRocha29/FamFlow.git
cd FamFlow
```

## Usage

To see the available commands:

```bash
python famflow.py --help
```

For help with a specific command:

```bash
python famflow.py <command> --help
```

## Project structure

A typical analysis may be organized into separate directories for the different stages:

```text
project/
├── sequences/
├── alignments/
├── models/
├── matrix/
├── ssn/
└── results/
```

The directory structure can be adapted according to the analysis.

## Design principles

FamFlow is based on three main principles:

**Reproducibility**
The same computational steps can be applied to different domains using a common workflow.

**Researcher control**
The tool assists with computational processing but does not determine the biological criteria used to construct a model.

**Integration**
Existing functionality, particularly from ROTIFER, is combined into a workflow rather than unnecessarily reimplemented.

## Applications

FamFlow can be used in analyses involving:

* Protein-family characterization
* Domain-level sequence analysis
* Custom HMM construction
* Remote homology detection
* Comparative genomics
* Sequence Similarity Network analysis
* Investigation of divergent protein families

## Project status

FamFlow is under development. The current implementation focuses on organizing sequence, alignment, and family-analysis workflows, with particular emphasis on applying the same computational procedures across multiple protein domains.

## Repository

[FamFlow — GitHub repository](https://github.com/MarcosRocha29/FamFlow?utm_source=chatgpt.com)

## Citation

If you use FamFlow in a publication, please cite the repository and the relevant software or methods used in the analysis.

---

**FamFlow**
*Workflow for protein family analysis and HMM model construction.*

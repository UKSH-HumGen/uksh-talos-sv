# experimental/

> **Status: Experimental — NOT integrated into the main workflow.**
>
> The code in this directory was prototyped during development of
> `uksh-talos-sv` but is **not invoked by the published pipeline**
> (`nextflow/talos.nf`) and was **not used to produce the results
> reported in the accompanying preprint**.
>
> It is preserved here for transparency and to give a starting point for
> anyone who wants to take these ideas further. Expect rough edges:
> module interfaces may not match the current `nextflow/` modules,
> inputs/outputs may be stale, and there is no test coverage.

## Contents

```
experimental/
├── nextflow/
│   └── modules/
│       └── talos/
│           ├── MergeStructuralVariants/   # alternative SV merging strategy
│           └── RescueCompoundHet/         # SNV+SV compound-het rescue
└── scripts/
    └── rescue_compound_het.py             # helper used by RescueCompoundHet
```

### `MergeStructuralVariants`
An alternative merge step that aggregates SVs across samples before the
main pipeline. Superseded in the active workflow by genotype-based
multi-sample SV assignment inside `RunHailFilteringSv` / downstream
steps. See `TALOS_SV_MODIFICATIONS.md` for rationale.

### `RescueCompoundHet`
A post-hoc step that pairs small variants and structural variants in the
same gene to flag potential SNV+SV compound heterozygous events.
Promising on paper but never validated end-to-end in this release.

## If you want to use this code

You are on your own. Wire the modules back into `nextflow/talos.nf`,
match their I/O against the current channels, and (please) add tests.

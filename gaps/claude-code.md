# Prompt for Code Review and Consistency Check

## Role

You are auditing the **GAPS** codebase (GPU-Amplified Parton Shower) for CPU↔GPU consistency. GAPS is a Monte Carlo event generator for parton showers and matching, with two parallel implementations that must produce bit-identical results when given the same seed:

- `gaps/cpu-shower/` — reference C++ implementation (CPU)
- `gaps/gpu-shower/` — CUDA C++ implementation (GPU)

**`cpu-shower` is always the source of truth for physics correctness.** When `gpu-shower` disagrees with `cpu-shower`, the GPU implementation is presumed wrong — unless the GPU version is clearly more correct (e.g. a bug fix, a more accurate constant, a corrected formula). Flag that explicitly as a "backend-ahead" issue: describe what the backend has that CPU lacks and why you believe the backend is right, and leave the decision to the reviewer. Do not silently let a backend correction go unrecorded.

All changes you make will be reviewed by the user before committing, so you can act with reasonable confidence — but you must clearly explain every change so the reviewer can judge it quickly.

## Goal

Audit the `matrix/`, `shower/`, `hadronisation/`, and `observables/` subdirectories in both trees to verify line-by-line consistency. This consistency is both a correctness requirement (identical numerical output) and a pedagogical one (new developers should be able to read the two side-by-side to learn GPU porting via CUDA).

Secondary goal: add docstrings to any function that lacks one, in either implementation.

Tertiary goal: flag pedagogy gaps — places where a new developer reading the code side-by-side would be confused or lose the physics thread.

**Do not refactor.** Changes must be minimal and targeted at consistency only. Do not rename symbols, reorganize files, modernize syntax, or clean up working code. If you notice something that could be improved but is not a consistency issue, ignore it.

## File pairing

**Drive the audit from `cpu-shower`.** For each source file in `cpu-shower/{matrix,shower,hadronisation,observables}/`, derive its GPU counterpart using this mapping:

| `cpu-shower` extension | `gpu-shower` extension |
| --- | --- |
| `.cpp` | `.cu` |
| `.h` | `.cuh` |

The directory structure mirrors exactly: a file at `cpu-shower/shower/src/shower.cpp` pairs with `gpu-shower/shower/src/shower.cu`.

**Pairing algorithm — follow this order for every file:**

1. Compute the expected GPU path from the CPU path using the table above.
2. Check whether the expected path exists.
   - If it exists: it is the counterpart. Proceed to the function-level audit.
   - If it does not exist: look for a file with the same base name and a variant suffix (e.g. `shower_kernel.cu`, `shower_gpu.cu`). If found, note the name deviation as a trivial observation but treat it as the counterpart. If nothing is found at all, record the CPU file as unpaired and continue.
3. Within a matched pair, enumerate functions and pair them. The pairing is not always 1:1 — see Non-1:1 pairings below.

**Non-1:1 pairings.** The GPU side commonly decomposes a single CPU function into a kernel/dispatch plus one or more device/inline helpers, or fuses several small CPU functions into a single kernel. When this happens:

- Treat the combined backend code (kernel + its helpers, or the fused kernel) as the counterpart of the corresponding CPU function(s).
- Check consistency across the union: every statement in the CPU function should appear somewhere in the combined backend code, and vice versa.
- Do not report these as missing or unpaired; they are paired, just with different decomposition.
- Note the decomposition pattern in the report (e.g. "CPU `evolve_shower` ↔ CUDA `evolve_shower_kernel` + `__device__ step_evolution`").

## Comment and docstring porting

This is the primary consistency task. Comments and docstrings are as load-bearing as code: a physics comment missing from the GPU backend is a consistency defect, not cosmetic drift.

**Porting directions — apply in this priority order:**

1. **CPU → GPU (mandatory).**  
   Every inline comment, section header, and docstring in `cpu-shower` must appear verbatim in `gpu-shower`. If the GPU backend is missing any of these, copy it. This is always a trivial fix.

2. **Backend → CPU (docstrings only, see Docstring rules).**  
   The only case where you write to `cpu-shower` is to add a docstring that exists in the GPU backend but not in CPU. This is still a trivial fix — it is a pure documentation addition with no logic change.

3. **Never port GPU-specific comments back to CPU.**  
   Backend-specific annotations (thread indexing, memory model notes, launch configuration rationale) stay in `gpu-shower`.

**What counts as a comment that must be ported:**

- Docstrings (any block comment above a function).
- Inline comments explaining a physics step, a formula, or a non-obvious value.
- Section divider lines and their labels.
- Any comment that explains *why* a line of code exists (intent, physics rationale, workaround).

**What does not need porting:**

- Comments explaining a GPU-specific adaptation that has no analogue in CPU (e.g. `// use __expf for device-side performance`).
- Anything explicitly flagged as backend-specific in the comment itself.

**Pedagogy exception.** If a comment exists in `cpu-shower` but was not ported to the GPU backend, it is simultaneously a consistency defect and a pedagogy gap. Fix it (port the comment) and also list it under Pedagogy gaps in the report so the reviewer understands the context.

## CPU → GPU function mapping

For every function in `cpu-shower`, there must exist a counterpart in `gpu-shower`.

**In `gpu-shower`**, the counterpart is one of: a `__device__` function, a `__host__ __device__` function, or a CUDA kernel (`__global__`) that wraps the logic.

The body of each counterpart must match the CPU version line-for-line wherever physically possible, including:

- function signatures (modulo backend-required qualifiers, pointer vs. reference for device memory),
- control flow and statement ordering,
- variable names,
- numerical constants and physics formulas,
- all physics comments and docstrings — verbatim.

**Legitimate divergences from CPU (do NOT flag):**

- CUDA qualifiers (`__device__`, `__host__`, `__global__`).
- Memory layout changes required by the backend (SoA vs AoS, raw pointers replacing STL containers).
- Launch configuration, thread indexing (`blockIdx`, `threadIdx`), and associated boilerplate.
- `std::` replaced by device-compatible equivalents (e.g. `sqrt`).
- Anything explicitly commented as a backend-specific adaptation.

## What to check for each function pair

1. **Signatures & logic** — same arguments (modulo backend type differences), same control flow, same operations in the same order.
2. **Comments & docstrings** — present in both implementations, textually identical (per the porting rules above).
3. **Variable names & statement ordering** — identical naming and sequencing.
4. **Numerical constants & formulas** — every literal, coefficient, and expression matches exactly.
5. **Pedagogy** — see Pedagogy check section.

Also flag: functions that exist in `cpu-shower` but are missing in `gpu-shower`, and functions in the GPU backend with no CPU counterpart (after accounting for decomposition patterns).

**Backend-ahead corrections.** If the GPU backend has a more correct constant, an additional guard, or a fixed formula that CPU lacks, do not silently treat it as a CPU-wins situation. Flag it in the Backend-ahead issues subsection.

## Docstring rules

If a function lacks a docstring in either implementation:

- **CPU has a docstring, GPU does not:** copy CPU docstring verbatim to the GPU side. Trivial fix.
- **GPU has a docstring, CPU does not:** copy GPU docstring verbatim to CPU. Trivial fix.
- **Neither has one:** write a new docstring and add it identically to both. Focus on *what* the function does and the *physics* it represents, not implementation mechanics. If you are unsure about the physics intent and would be guessing, skip it and list the function under "docstrings skipped" in the report.

Match the project's existing docstring style — read a few existing ones before writing any new ones.

After this pass, every paired function must have the same docstring across both implementations.

## Pedagogy check

This codebase is pedagogical: new developers read `cpu-shower` and `gpu-shower` side-by-side to learn how a physics kernel is ported to CUDA.

**Flag as a pedagogy gap if:**

- A non-obvious physics step has no comment explaining *why* it is done (Sudakov veto weight, kinematic bound, splitting-function normalization).
- A backend-specific transformation (pointer arithmetic, SoA index, CUDA warp operation) is not annotated with what CPU concept it replaces.
- A docstring describes *what* the function does but not *what physics it represents*.
- A variable meaningful in the physics context is used without any comment at first use (`z`, `t`, `kt2`, etc.).
- A magic number is present with no comment. This includes physics constants expressed as plain float literals.
- A comment exists in `cpu-shower` explaining a subtlety but was not ported to the GPU backend — a developer reading only the GPU version would miss it. (Also fix this under the comment porting rules above.)

**Do not flag:** boilerplate obvious to any C++/CUDA developer; comments that are already sufficient; conventionally terse physics variable names (`Q2`, `x`, `mu`).

**Do not add or rewrite comments during the audit pass.** Pedagogy gaps are report-only. The exception: if a pedagogy gap is also a missing-comment consistency issue (the comment exists in CPU but was not ported), fix it as a trivial fix and note it in the pedagogy section.

## Fix vs. report policy

Since the user reviews everything before committing, err toward making the fix and explaining it clearly rather than deferring.

**Fix directly and note briefly:**

- Missing or differing physics comments and docstrings — copy verbatim.
- Adding new docstrings per the rules above.
- Typos, renamed local variables, reordered independent statements.
- Whitespace, formatting, and comment-only drift.
- Missing `const`, obviously intended constant values where the CPU value is clearly correct.
- README additions for missing CLI arguments (append only, same style).

**Fix with detailed explanation, or report without fixing if uncertain:**

- Any change to physics logic, formulas, or numerical constants used in a computation.
- Differing control flow (branches, loop bounds, early returns).
- Missing functions or kernels — if the port is mechanical, write it; if it requires design decisions (launch config, memory layout), report instead.
- Signature mismatches beyond backend-qualifier adaptation.
- Any `interface/` or `rungaps` inconsistency that touches argument parsing, dispatch logic, or public API surface.

For complex fixes, state: what changed, why the reference version is correct, and what numerical impact (if any) you expect. If you cannot confidently state the numerical impact, report rather than fix.

**Default to reporting rather than fixing when any of these apply:**

- The divergence involves floating-point specifics: `float` vs `double`, fused multiply-add, reciprocal or rsqrt approximations, `-ffast-math` behavior, or operation ordering in a sum/product.
- The GPU version uses a CUDA intrinsic (`__expf`, `__logf`, `__sinf`) where CPU uses the standard library — this is usually intentional.
- The divergence has been stable across multiple commits (`git log -p` or `git blame`) — longevity is weak evidence of intentionality.
- There is any comment near the divergent code that might be justifying the difference.
- The function involves random number generation, state indexing into an RNG stream, or anything order-dependent across threads.

In these cases, report with a note that it *might* be intentional. A false fix here silently breaks the bit-identical guarantee.

## Scope: additional files

Beyond the two subdirectories, audit:

- `gaps/interface/` — check for internal consistency across all files, and verify that the interface is consistent with how `cpu-shower` and `gpu-shower` consume or expose it.
- `rungaps` — check for consistency with `interface/` and the two shower implementations.
- `README.md` — cross-reference the documented CLI arguments against those actually accepted by `rungaps` and `interface/`. Add any missing arguments. Do not remove or rewrite existing entries; only append, in the same style.

Shared headers outside the two subdirectories are **out of scope** for physics-logic changes. Report shared-header issues rather than fixing them, noting the header path so the user can address them separately.

## Workflow

Process in this order: `matrix/`, `shower/`, `hadronisation/`, `observables/`, then `interface/` + `rungaps` + `README.md` as a single final section. For each subdirectory:

1. **Orchestrator (you) builds the pairings.** List all files in both trees and build the CPU↔CUDA file pairings using the Pairing algorithm above (including any kernel + helper decompositions). This step stays on the orchestrator — do not delegate pairing, because it needs a whole-subdirectory view to spot variant-suffixed and unpaired files.
2. **Dispatch one subagent per CPU↔GPU file pair** (see Parallel subagent execution below). Each subagent owns exactly one pair, performs the full function-level audit on it, applies its own trivial and confident complex fixes, and returns a structured report fragment.
3. **Orchestrator aggregates.** Collect every subagent's fragment, merge them into a single subdirectory section using the report template, and append that section to `AUDIT_REPORT.md` before moving on to the next subdirectory.

Unpaired CPU files (no counterpart found in step 1) and shared-header issues do not get a subagent — record them directly in the aggregated section.

The `interface/` + `rungaps` + `README.md` scope is cross-cutting rather than pair-structured; audit it yourself on the orchestrator rather than delegating per pair.

## Parallel subagent execution

The per-pair audit parallelizes cleanly: each CPU↔GPU pair touches a disjoint set of files, so subagents editing different pairs can never write-conflict on source files. Use the `Agent` tool (`general-purpose` subagent) to run pairs concurrently — dispatch a batch of pairs in a single response so they execute in parallel.

**Each subagent's task prompt must include:**

- The exact CPU path and GPU path for its one pair, plus any kernel/helper decomposition already identified for it.
- A pointer to this instruction file (`gaps/claude-code.md`) as the governing rules, with the load-bearing constraints restated inline so the subagent does not have to infer them: `cpu-shower` is the source of truth; the only permitted CPU edits are docstring additions/copies; do not refactor or rename; default to reporting (not fixing) on floating-point, CUDA-intrinsic, and RNG-order divergences; port comments and docstrings verbatim; any new code the subagent writes must be snake_case and well commented, and existing non-snake_case symbols are reported (category `naming-convention`), never renamed.
- The instruction to apply all checks from "What to check for each function pair," make trivial and confident complex fixes directly to its two files, and accumulate uncertain items for its fragment.

**Each subagent must NOT:**

- Touch any file outside its assigned pair. This keeps parallel edits conflict-free and is what makes the fan-out safe.
- Write to `AUDIT_REPORT.md`. Report writes are serialized on the orchestrator to avoid concurrent-write races; the subagent returns its findings as its final message instead.

**Each subagent returns a fragment** covering only its pair, using the per-pair rows of the report template: trivial fixes applied, complex fixes applied (with what/why/expected numerical impact), issues reported, backend-ahead issues, docstrings added/skipped, pedagogy gaps, non-1:1 pairings noted, and unpaired functions within the pair. The orchestrator does not re-audit — it trusts each fragment and merges them, summing the counts into the subdirectory Summary table.

**Concurrency guidance.** Dispatch pairs within a subdirectory in parallel, but keep subdirectories sequential (append each subdirectory's aggregated section to `AUDIT_REPORT.md` before starting the next) so the report is written incrementally and progress is never lost. If a subagent fails or returns an incomplete fragment, re-dispatch that single pair rather than restarting the subdirectory.

## Output: write to `AUDIT_REPORT.md`

Create `AUDIT_REPORT.md` at the repository root. Write to it incrementally after each section so progress is not lost. Start with a heading and a run metadata block (date, commit hash, subdirectories audited). Then append one section per subdirectory using this template:

````markdown
## `<subdirectory>/`

**Summary**

| Metric | Count |
|---|---|
| Files audited (CPU↔CUDA pairs) | N |
| Trivial fixes applied | M |
| Complex fixes applied | P |
| Issues reported (not fixed) | K |
| Backend-ahead issues flagged | B |
| Docstrings added | D |
| Docstrings skipped | S |
| Pedagogy gaps flagged | G |

### Trivial fixes applied

**`<cpu_path>` ↔ `<gpu_path>`**
- <one-line description of fix>

### Complex fixes applied

**`<cpu_path>` ↔ `<gpu_path>`**
- **`<function_name>`**
  - *What changed:* <short description>
  - *Why:* <reference to the correct version + reasoning>
  - *Expected numerical impact:* <none | describe>

### Issues reported (require your review)

**`<cpu_path>` ↔ `<gpu_path>`**
- **`<function_name>`** — *category: <missing | signature | logic | constant | control_flow | numerics | other>*
  - CPU (line X): <short excerpt or description>
  - GPU (line Y): <short excerpt or description>
  - Why I didn't fix it: <one line>

### Backend-ahead issues (GPU may be more correct than CPU)

**`<gpu_path>`**
- **`<function_name>`** — *category: <constant | formula | guard | other>*
  - GPU (line X): <what the backend has>
  - CPU (line Y): <what CPU has, or "missing">
  - Why this may be a backend fix: <one line>

### Docstrings added

- Copied CPU → CUDA: <count>
- Copied CUDA → CPU: <count>
- Written new (both sides): `<function_name>`, ...
- Skipped (physics intent unclear): `<function_name>`, ...

### Pedagogy gaps

> Locations where a new developer reading the code would lose the physics thread. No fixes applied — flagged for the reviewer.

- **`<file_path>` line X** — *category: <missing-physics-comment | missing-port-annotation | docstring-what-not-why | unexplained-variable | magic-number | comment-not-ported | naming-convention>*
  - <one-line description of what is missing and why a newcomer would be confused>

### Non-1:1 pairings noted

**CPU↔CUDA**
- CPU `<function>` ↔ CUDA `<kernel>` + `<helpers>`

### Unpaired functions

**CPU↔CUDA**
- CPU-only: `<function_name>`, ...
- CUDA-only: `<function_name>`, ...

### Shared-header issues (out of scope, reported for reference)

- `<header_path>`: <description>
````

For the `interface/` + `rungaps` + `README.md` section, use the same template but adapt the "Files audited" rows to list interface files, and add:

````markdown
### README.md — missing arguments added

- `<argument>`: <one-line description, matching existing README style>

### README.md — arguments present in README but not found in rungaps/interface (possible stale docs)

- `<argument>`: <note>
````

After all sections, append a final `## Overall summary` with totals and a short prose summary highlighting: the most important issues for the reviewer to look at first; patterns noticed across subdirectories; divergences flagged as possibly intentional; backend-ahead issues where the GPU disagrees with CPU (most likely genuine CPU bugs); and highest-priority pedagogy gaps (recurring missing physics explanations, or backend-specific transformations never annotated anywhere).

## Code style for anything you write

Whenever the audit produces new code — a missing-port kernel or function, a new local variable, a new docstring or comment — it must follow the project convention:

- **snake_case** for all identifiers: functions, kernels, variables, and parameters. CUDA qualifiers and built-ins (`__global__`, `blockIdx`, `threadIdx`) are exempt, as are physics symbols that are conventionally written otherwise (`Q2`, `kt2`). When in doubt, match the casing of the surrounding existing code.
- **Well commented:** every new function or kernel gets a docstring (per the Docstring rules), and any non-obvious physics step, formula, or magic number inside it gets an inline comment explaining *why*. New code is held to the same pedagogy bar the audit flags elsewhere — do not introduce a pedagogy gap while fixing one.

This applies to **new** code only. It does **not** license renaming existing symbols — the no-refactor / no-rename constraint below always wins. If an existing symbol violates snake_case, leave it untouched and record it as a report-only pedagogy observation (category `naming-convention`); never rename it, and never let a CPU↔GPU casing mismatch drive you to rename a paired symbol — flag it instead.

## Constraints

- All code is snake_case and well commented, per Code style for anything you write above. This never overrides the no-rename rule: existing symbols are reported, not renamed.
- The only permitted modifications to `cpu-shower/` are: copying an existing GPU docstring over, or adding a new docstring to a previously undocumented function. No logic changes to CPU code, ever.
- Do not modify any files outside the two audited subdirectories and the `interface/` + `rungaps` + `README.md` scope — report shared-header issues instead.
- Do not change GPU code in ways that alter numerical behavior unless you can clearly justify the change with reference to the correct source, and be especially cautious around floating-point and CUDA intrinsics.
- Preserve existing backend-specific adaptations in `gpu-shower`; don't fix them toward the CPU version unless CPU is clearly right.
- Match the project's existing comment and docstring style rather than imposing a new one.
- Do not refactor, rename, or reorganize code. Keep the diff small.

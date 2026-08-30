# 11 · Making the analysis reproducible

Reproducibility here is not a slogan; it is four concrete design decisions.

**1 · Everything runs from a workflow, not from a shell history.** `workflow/Snakefile` declares the dependency
graph; Snakemake reruns only what changed and refuses to continue past a failed step. `snakemake -n` prints the
plan without executing it. The same workflow file describes both execution modes.

**2 · Software versions are pinned by environment files.** The `environment/*.yaml` files are resolved by conda
from Bioconda, so the same tool versions can be installed on any machine. `scripts/check_tools.py` verifies that
every required binary is present *and executable* — Bioconda occasionally resolves a Linux build into a macOS
environment, which otherwise fails deep inside a run with an opaque exit code 126.

**3 · The public example is real, small and self-contained.** `scripts/fetch_example_reference.py` fetches ~230 kB
of genuine GRCh38 sequence around six driver-gene hotspots; `scripts/simulate_amp.py` builds an anchored-multiplex
library over it, spiking known mutations (KRAS G12C, BRAF V600E, TP53 R175H, EGFR L858R, PIK3CA H1047R,
PTEN R130*) at defined allele fractions, with PCR duplicates and sequencing errors, plus a truth file. Because
the coordinates are real, live annotation works on the example exactly as on a patient sample — and because the
reference is tiny, the whole pipeline runs on a laptop in under a minute. `tests/test_pipeline_outputs.py` then
asserts that the spiked variants are recovered, classified somatic, annotated with the right protein change and
tiered correctly. That is the difference between "the workflow ran" and "the workflow is right".

**4 · Heavy steps are reproducible too.** Alignment against the human genome runs on usegalaxy.eu; Galaxy
workflows export as a `.ga` file that anyone can import and rerun, and the tool versions are recorded in it. The
exported workflows and the exact parameters live in `galaxy/`.

**Provenance.** Every JSON output carries a `_provenance` block with the UTC timestamp, the git commit, whether
the working tree was dirty, the Python version and the command line. If a number in the report is questioned,
the run that produced it can be identified.

**Why not Nextflow/nf-core.** `nf-core/sarek` and `nf-core/oncoanalyser` are the production-grade equivalents and
are the right answer in a laboratory with a cluster: oncoanalyser's own documentation allocates 72 GB of RAM and
375–750 GB of disk to its standard steps, and its targeted mode needs a panel reference built from ≥20 samples.
Neither fits an 8 GB laptop or a single sample, which is why this project uses Snakemake plus Galaxy — and says
so, rather than pretending the constraint does not exist.

# Verifying parallelism + skip behaviour

Use these recipes on the Linux server while the pipeline is running to
confirm that the output-aware skip and (eventually) parallel workers are
actually engaging.

---

## Output-aware skip (PR 1)

### Quick check: is the skip path firing?

After a successful step-10 (CCI) run, re-launch the same command and
look for `↻ X skipped: ... already exists` lines in the log:

```bash
./Run_Scripts/Run_Giotto_Pipeline.sh \
  --sample-steps 10_cci --samples CART_T0_S1 \
  2>&1 | tee /tmp/skip_check.log

grep -E "↻|skipped: outputs already exist|skipped: checkpoint exists" /tmp/skip_check.log
```

Expected output:
```
  ↻ InSituCor section skipped: CART_T0_S1_insitucor_module_stats.csv already exists. ...
  ↻ LIANA section skipped: CART_T0_S1_liana_aggregate.csv already exists. ...
  ↻ NicheNet section skipped: CART_T0_S1_nichenet_comparison_summary.csv already exists. ...
  ↻ nnSVG section skipped: CART_T0_S1_nnSVG_results.csv already exists. ...
```

### Force regeneration

```bash
./Run_Scripts/Run_Giotto_Pipeline.sh \
  --sample-steps 10_cci --samples CART_T0_S1 --overwrite
```

The banner should print:
```
⚠ --overwrite active: existing checkpoints / section outputs will be regenerated rather than skipped.
```

And no `↻ ... skipped` lines should appear.

### Verify a sentinel file is actually checked

```bash
# Move a sentinel out of the way
mv Output/Sample_CART_T0_S1/10_CCI_Analysis/svg/CART_T0_S1_nnSVG_results.csv /tmp/

# Re-run — only nnSVG should re-execute; the others stay skipped
./Run_Scripts/Run_Giotto_Pipeline.sh --sample-steps 10_cci --samples CART_T0_S1
```

---

## smiDE thread usage (PR 1: `smi_de(nCores=)` only)

`smide_ncores: 4` in `config.yaml` reaches `smiDE::smi_de(nCores=)` at
the GLM-fit step. Other smiDE calls (`pre_de`, `overlap_ratio_metric`,
`results`) are single-threaded internally regardless of this value.

Confirm at runtime:

```bash
# Identify the parent R process
PID=$(pgrep -f "Rscript.*CosMx_pipeline" | head -1)

# Live thread view — when smi_de() is fitting, expect ~smide_ncores busy threads
top -H -p ${PID}

# Snapshot active threads
ps -L -p ${PID} -o pid,tid,pcpu,pmem,comm | head -20

# Confirm BLAS / OpenMP env (BLAS threads can compete with smi_de's nCores)
apptainer exec ${SIF} Rscript -e \
  "cat('OMP_NUM_THREADS=', Sys.getenv('OMP_NUM_THREADS'),
       '\nOPENBLAS_NUM_THREADS=', Sys.getenv('OPENBLAS_NUM_THREADS'), '\n', sep='')"
```

Note: between cell-type iterations (`pre_de`, `overlap_ratio_metric`,
`results`) only one thread is active. The 4-thread spike happens during
the actual GLM fit inside `smi_de()`.

---

## Coming in PR 2

- **Sample-level parallelism** (`--workers N`) — multiple samples run
  concurrently in separate R sessions via `future.apply::future_lapply`.
- **smiDE per-cell-type parallelism** (`smide_inner_workers: N`) — each
  cell type's `pre_de + smi_de + results` triple is dispatched to its
  own future worker. Compounds with sample workers (`N_outer × N_inner`
  total processes — pick both based on `nproc; free -g`).

For PR 2, additional verification recipes will cover:

```bash
# Snapshot active workers
ps --ppid ${PID} -o pid,pcpu,pmem,rss,cmd | grep -v grep
# Expect N child R processes when --workers N or smide_inner_workers: N
```

---

## Server pre-flight (always)

Before any multi-worker run, check what the server is actually doing:

```bash
nproc                # total logical cores
free -g              # available RAM (look at "available" column)
uptime               # 1/5/15 min load averages
who                  # other users active
```

Pick `--workers` and `smide_inner_workers` such that
`workers × smide_inner_workers + 2` is well under available cores AND
`workers × ~80 GB` fits in available RAM.

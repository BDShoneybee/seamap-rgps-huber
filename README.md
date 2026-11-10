# Seamap RGPS Huber trend study

Research reconstruction and reproducibility package for **Causal robust trend preprocessing of marine relative-GPS range and bearing: reproducible evidence and operational limits**.

Authors: Menghui Fan, Shouchuan Fang and Fujiang Liu.

## Reproduce

From the repository root:

```sh
python -m pip install -r requirements.txt
python paper/scripts/experiments.py
python paper/scripts/extended_experiments.py
python paper/scripts/archive_diagnostics.py
python paper/scripts/make_tables.py
python paper/scripts/verify_native.py
cd paper
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

The native verification additionally requires a C++14 g++ compiler. Python >=3.10 is expected; the retained outcomes were produced with Python 3.14 and NumPy as recorded in `paper/results/verification.json`. Matplotlib generates vector PDF and 300 dpi PNG figures. A current pdfLaTeX installation with amsmath, amssymb, graphicx, hyperref and lineno compiles the manuscript. Line numbering is available in normal TeX distributions. No Python PDF-generation substitute is used: equations and manuscript are authored in LaTeX.

`experiments.py` uses the original source tree if present beside `paper/`. Otherwise it consumes the extracted public NPZ archives, preserving the same input values, gaps and elapsed time. Existing extracted log checks and event counts support the log audit without redistributing mixed application logs. The default script does not infer field truth or equipment precision.

## Inputs and definitions

- `paper/data/range_archive.npz`, `bearing_archive.npz`: each contains arrays `t` (elapsed seconds), `y` (input range in metres / bearing in radians, missing as NaN) and `archived_prediction` (historical output, **not truth**). Each channel has 8830 epochs and 6741 observed values.
- `paper/results/*_replay.csv`: human-readable inputs and rerun outputs. Columns `time`, `input`, `trend`, `huber_kf`. Units match the corresponding archive.
- `paper/results/monte_carlo_trials.csv`: independent scalar run-level RMSE, MAE and P95 absolute error; normalized units; 60 seeds per scenario.
- `simulation_summary.csv`: means and 95% run-level intervals (`ci95` is half width).
- `paired_bootstrap.csv`: Trend minus comparator RMSE, with paired run-level bootstrap interval.
- `injection_trials.csv`: rerun-to-rerun output divergence after injected errors; not truth error.
- `synthetic_position_trials.csv`: horizontal RMSE in metres against **simulated** truth.
- `step_response_trials.csv`: post-step RMSE and first entry into the 10% band, not sustained settling time.
- `gap_subset_trials.csv`: separate errors for observed, short-gap and long-gap epochs in the combined scenario.
- `detection_trials.csv`: detection counts against injected synthetic outlier labels only.
- `logged_equation_checks.csv`: arithmetic discrepancies obtained from rounded logged priors and logged accepted-update parameters.
- `native_parity.json`: portable C++ / Python agreement; this is not full production-application equivalence.

The input source is ship-collected Seamap RGPS equipment, as confirmed by the authors. Some diagnostic folders in the retained project refer to indoor tests; live processing versus later replay is not independently documented. Device generation, site, sea state, receiver precision and independent reference trajectory are unavailable. Range absolute hour and bearing absolute hour are not jointly established in the extracted archive. Channel replay uses independently recovered elapsed times. Available observations within the first 300 archive epochs determine normalization; evaluation begins at archive index 300. Parameters are fixed, not selected on evaluation outcomes.

## Algorithm scope

The primary algorithm reconstructs the scalar Huber core in DLPPreProc, with an explicit six-observation warmup anchored at the latest initialization epoch. The Qt/Eigen polynomial wrapper, configured equipment-accuracy floor and application reset paths are not claimed to be reproduced. Missing epochs query the latest state: extrapolate up to 8 s, freeze until 100 s, unavailable after 100 s. Missing queries do not mutate state. The portable native code applies increasing timestamp validation. Angular ring residuals in Python are an explicit optional extension, not silently attributed to the scalar production code.

`No-alternation` and `No-gate` are matched component ablations. `Huber-KF` is the precise single-pass weighted linear KF documented in the manuscript, not a reproduction of a named published HCKF. Raw/hold and median-7 use only available historical samples. The KF process intensity is 0.0025 and normalized measurement variance is 1. No claim is made that all comparators have globally optimal tuning.

## Interpretation

The trend method improves slow-trend synthetic results with isolated outliers but is not consistently best. Ship bearing becomes smoother while next-observation disagreement increases. No independent reference exists, so there is **no claim of measured absolute position/velocity accuracy**. The public record is a noisy experimental input, not a surrogate ground truth. See the manuscript for limitations and precise protocols.

## Manuscript status

The manuscript is an author-review draft. Corresponding author: Shouchuan Fang (fangshouchuan@cnpc.com.cn). No competing interests and the listed contributions were confirmed by the authors. Funding is left blank as requested. Final manuscript approval remains the responsibility of the authors. Author approval and submission are not asserted by this package. GitHub access provides reproducibility; it does not supply a dataset DOI or peer review.

## Rights

No new third-party software or raw-log license is invented. Publication permission was requested by the data provider in this work session. No blanket software/data license has been selected on behalf of the authors. Reuse beyond inspection and reproduction should be agreed with the authors. `lineno.sty` is third-party LaTeX support distributed under its own license, not original research software.

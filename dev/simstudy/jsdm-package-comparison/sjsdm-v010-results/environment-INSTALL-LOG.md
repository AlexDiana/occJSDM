# sjSDM v0.1.0 isolated environment: install log

Installed on 13 October 2026 by Claude Code at Doug's request, for the TODO item "Check that sjSDM v0.1.0 gives the same fit as the recorded one". Lesson 5 tells readers to install sjSDM from Doug's fork release v0.1.0, but the recorded sjSDM fits ran under fork v0.2.1 (commit d2ca508) with `SJSDM_MOJO_BACKEND=0`, that is, PyTorch on the CPU. This directory only prepares the environment for rerunning `sjsdm-multistart.R` (seed 26092811) under v0.1.0. Installation and a smoke test were done; the comparison itself has not been run.

Everything lives in this directory. Nothing was installed into the R 4.6 user library or the system R library, no conda environment was created, `install_sjSDM()` was not called, and the fork checkout was not modified.

## Layout

- `source/`: the `git archive` export of the fork's `sjSDM/` subdirectory at tag v0.1.0.
- `venv/`: the Python virtual environment.
- `library/`: the dedicated R 4.6 library (reticulate, sjSDM, and the sjSDM dependencies missing from the R 4.6 user library).
- `run-r`: the launcher. Use it for every R process that should see this environment.
- `smoke-test.R` and `smoke-test.log`: the smoke test and its output.
- `pip-freeze.txt` and `r-session-info.txt`: version records.

## Fork source

- Fork checkout: /Users/douglasyu/src/s-jSDM (branch master, HEAD 7441c30, clean `git status` before and after this install).
- Tag v0.1.0 resolves to commit e7ee5826a1a1fbef951f140c4feb186b22f84e11.
- The R package inside reports sjSDM 1.0.7, the same R package version as fork v0.2.1. The fork release number and the R package version are separate, so record both.
- The bundled Python module `sjSDM_py` has version 0.1.3 in its setup.py, which declares only numpy. The tag's `_environment.yml` lists pyro-ppl, madgrad, `torch_optimizer`, torch and pytest, all unpinned.

Export command:

``` sh
git -C ~/src/s-jSDM rev-parse 'v0.1.0^{commit}'
git -C ~/src/s-jSDM archive v0.1.0 sjSDM | tar -x -C source
```

## Python

- Base interpreter: /Users/douglasyu/.local/bin/python3.12, CPython 3.12.13, the uv-managed python-build-standalone build (Clang 22.1.3), base prefix /Users/douglasyu/.local/share/uv/python/cpython-3.12.13-macos-aarch64-none.
- Venv executable: `venv/bin/python` in this directory. The venv has no pip; use `uv pip` against it.
- Package source: PyPI only.

Commands, run from this directory:

``` sh
~/.local/bin/uv venv --python /Users/douglasyu/.local/bin/python3.12 --no-project venv
~/.local/bin/uv pip install --python venv/bin/python --index-url https://pypi.org/simple \
  "torch==2.5.1" "numpy==2.5.2" "pyro-ppl==1.9.1" "pyro-api==0.1.2" "opt-einsum==3.4.0" \
  "tqdm==4.70.0" "madgrad==1.3" "torch-optimizer==0.3.0" "pytorch-ranger==0.1.1" \
  "pytest==9.1.1" "pluggy==1.6.0" "iniconfig==2.3.0" "packaging==26.3" "pygments==2.21.0" \
  "filelock==3.32.3" "networkx==3.6.1" "typing-extensions==4.16.0" "jinja2==3.1.6" \
  "markupsafe==3.0.3"
```

Installed packages (the full list, identical to `pip-freeze.txt`):

- torch 2.5.1
- numpy 2.5.2
- pyro-ppl 1.9.1
- pyro-api 0.1.2
- madgrad 1.3
- torch-optimizer 0.3.0
- pytorch-ranger 0.1.1
- tqdm 4.70.0
- opt-einsum 3.4.0
- pytest 9.1.1
- pluggy 1.6.0
- iniconfig 2.3.0
- packaging 26.3
- pygments 2.21.0
- filelock 3.32.3
- networkx 3.6.1
- typing-extensions 4.16.0
- jinja2 3.1.6
- markupsafe 3.0.3
- sympy 1.13.1
- mpmath 1.3.0
- fsspec 2026.9.0
- setuptools 84.0.0

### Match to the original environment

The Python package versions match the environment the recorded v0.2.1 fits used, except where noted below. They were not chosen fresh.

Evidence. DESIGN.md (line 27) says the v0.2.1 runs used the existing working Python environment, Python 3.12.13 and PyTorch 2.5.1, reached through the repository's reticulate venv. No pip freeze or `py_config()` record of the Python packages was saved with the results: the fit records store only the Python path, and their saved `sessionInfo()` covers R only. The fork's AGENTS.md identifies that venv as `work/mojo-backend/reticulate-venv` in the fork checkout, a `--system-site-packages` wrapper over the pixi environment `work/mojo-backend/.pixi/envs/default`. Contrary to the brief, that environment still exists and still imports. Its `pixi.lock`, its site-packages metadata and a read-only import check (`python -I -B`, so no bytecode was written into the checkout) give torch 2.5.1, numpy 2.5.2, pyro-ppl 1.9.1, pyro-api 0.1.2, madgrad 1.3, torch-optimizer 0.3.0, pytorch-ranger 0.1.1, tqdm 4.70.0, opt-einsum 3.4.0 and pytest 9.1.1. Its package files are dated 22 or 28 August 2026, before the September runs, so it is very likely the exact set those fits used. All of these versions were pinned here, along with the matching versions of the small helper packages (filelock, networkx, typing-extensions, jinja2, markupsafe, packaging, pygments, pluggy, iniconfig).

Not matched:

- torch 2.5.1 is the PyPI wheel. The original was the conda build `pytorch-2.5.1-py3.12_0` from the pytorch channel. Their `torch.__config__.show()` output was compared: same BLAS and LAPACK (Accelerate), same compiler flags, OpenMP parallel backend, MKL off, and the same default thread counts. The only differences are the compiler's path name and the spelling of `USE_CUDA` (0 versus OFF).
- The Python interpreter is the python-build-standalone 3.12.13 (Clang 22.1.3). The original was the conda-forge 3.12.13 (Clang 19.1.7).
- numpy 2.5.2 is the PyPI wheel, which links Accelerate. The original was the conda-forge build on the reference BLAS interface backed by OpenBLAS 0.3.34. sjSDM does its arithmetic in torch. It uses numpy for data transfer and for the random initial values of sigma, which come from numpy's random generator, not from BLAS.
- sympy 1.13.1 and mpmath 1.3.0 (original 1.14.0 and 1.4.1) are forced by the torch 2.5.1 wheel, which pins `sympy==1.13.1`.
- fsspec and setuptools are extra, because the torch wheel requires them.
- The original environment also contained Mojo 1.0.0, MAX 26.5.0 and unrelated tooling. These were not installed, since v0.1.0 has no Mojo code.

## R

- R 4.6.1 (2026-06-24), /usr/local/bin/Rscript, a symlink into /Library/Frameworks/R.framework. Doug decided on 13 October 2026 that new work uses R 4.6. The recorded v0.2.1 fits ran under R 4.5.0 with reticulate 1.46.0, according to the `sessionInfo()` saved in their fit records.
- Library path under `run-r`: this directory's `library/`, then /Users/douglasyu/Library/R/arm64/4.6/library, then the system library.
- Package source: CRAN cloud mirror (`cloud.r-project.org`) macOS arm64 binaries for the dependencies. sjSDM was installed from the exported source.

The recursive Depends, Imports and LinkingTo of sjSDM's Imports came to 67 non-base packages. Twelve were missing from the user and system libraries and were installed into `library/` with `dependencies = FALSE`:

- reticulate 1.47.0
- Metrics 0.1.4
- mathjaxr 2.0-0
- beeswarm 0.4.0
- qgam 2.0.0
- viridis 0.6.5
- RcppTOML 0.2.3
- here 1.0.2
- png 0.1-9
- doParallel 1.0.17
- foreach 1.5.2
- iterators 1.0.14

Then sjSDM itself (R package 1.0.7) went into the same library. Both load tests in `R CMD INSTALL` passed:

``` sh
lib=/Users/douglasyu/src/occJSDM/dev/simstudy/results/sjsdm-v0.1.0-env-20261013/library
Rscript --vanilla -e 'install.packages(c("reticulate","Metrics","mathjaxr","beeswarm","qgam","viridis","RcppTOML","here","png","doParallel","foreach","iterators"), lib = "'"$lib"'", repos = c(CRAN = "https://cloud.r-project.org"), dependencies = FALSE)'
env R_LIBS="$lib:/Users/douglasyu/Library/R/arm64/4.6/library" \
  RETICULATE_PYTHON=/Users/douglasyu/src/occJSDM/dev/simstudy/results/sjsdm-v0.1.0-env-20261013/venv/bin/python \
  KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /usr/local/bin/R CMD INSTALL --no-multiarch --library="$lib" source/sjSDM
```

All the other dependencies come from the R 4.6 user library. `r-session-info.txt` lists every dependency with the library it resolved from. No file in the user library has a modification date of 13 October 2026.

## Launcher

`run-r` exports the following and then runs /usr/local/bin/Rscript with the given arguments. Use `./run-r script.R` or `./run-r -e 'expr'`.

- `R_LIBS`: this directory's `library/`, then the R 4.6 user library.
- `RETICULATE_PYTHON`: this directory's `venv/bin/python`. v0.1.0's `.onLoad` only looks for a conda environment named r-sjsdm, inside a silent `try()`. There is no conda on this Mac, so that lookup fails quietly and `RETICULATE_PYTHON` decides. `py_config()` confirms this: "Python version was forced by RETICULATE_PYTHON".
- `KMP_DUPLICATE_LIB_OK=TRUE`.
- `SJSDM_MOJO_BACKEND=0`. v0.1.0 ignores it. It is set to match the earlier runs and because `sjsdm-multistart.R` checks for it.
- `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`. Under the launcher, torch reports 1 thread.

## Smoke test

All three steps passed.

Step 1, the Python imports:

``` sh
venv/bin/python -c 'import torch, pyro, madgrad, torch_optimizer, numpy; ...'
```

This printed torch 2.5.1, pyro-ppl 1.9.1, madgrad 1.3, torch-optimizer 0.3.0, numpy 2.5.2, pytorch-ranger 0.1.1, tqdm 4.70.0 and opt-einsum 3.4.0. A float64 matrix product worked, and MPS was reported available (it was not used).

Steps 2 and 3 are in `smoke-test.R`, with the full output in `smoke-test.log`:

``` sh
./run-r smoke-test.R
```

Step 2, loading sjSDM in R:

- `library(sjSDM)` attached sjSDM 1.0.7 from `library/sjSDM`, with reticulate 1.47.0 from `library/reticulate`.
- `reticulate::py_config()$python` is this directory's `venv/bin/python`, with libpython from the uv 3.12.13 install, numpy 2.5.2 and torch from the venv.
- The bundled module `sjSDM_py` imported from `library/sjSDM/python/sjSDM_py/__init__.py` and has `Model_sjSDM`.
- `check_installation()` found torch, `torch_optimizer`, pyro and madgrad. torch reports version 2.5.1 and 1 thread.

Step 3, two tiny CPU fits on data from `simulate_SDM(env = 3L, species = 5L, sites = 100L, seed = 42L)`:

- Fit A was the requested default fit: `env = linear(data = com$env_weights, formula = ~ .)`, `se = FALSE`, `iter = 20L`, `device = "cpu"`, `sampling = 100L`, `seed = 42L`. It returned a finite 5 by 4 coefficient matrix. Device cpu, dtype `torch.float32`, 20 history entries, final loss 2.415.
- Fit B used the argument set from `sjsdm-multistart.R`, but with 20 epochs and 100 draws: `binomial("logit")`, `dtype = "float64"`, `bioticStruct(df = 2, lambda = 0)`, `linear(lambda = 0)`, `step_size = 100L`, `learning_rate = 0.002`, and `sjSDMControl(optimizer = RMSprop(weight_decay = 0.0001), scheduler = 0, early_stopping_training = 0)`. v0.1.0 accepted every argument. It returned finite coefficients. Device cpu, dtype `torch.float64`, df 2, and `force_r(model$env_weights)` and `force_r(model$get_sigma)` have the expected shapes.

## Before running the v0.1.0 comparison

- The unmodified `sjsdm-multistart.R` will stop under v0.1.0. It asserts `settings$backend == "torch_tmp"` in R and `m._loss_function.__name__ == 'torch_tmp'` in its embedded Python. In v0.1.0 the same eager PyTorch logit loss is named `tmp`: fit B reports `_loss_function` name `tmp`. v0.2.1 renamed it `torch_tmp` when it added the Mojo switch. Run an adapted, immutable copy that accepts `tmp`. Apart from the name, the only change to that loss between the tags is the `int()` casts on the random-draw dimensions added in v0.2.1. Neither version applies TorchScript to it.
- The script also requires `start` to be in 4:12 (seed 26092811 is start 11) and `SJSDM_MOJO_BACKEND` to be 0. `run-r` sets the latter.
- Each Python start prints a SyntaxWarning from a pyro 1.9.1 docstring. It is harmless. It repeated here on every run because the installing shell had `PYTHONDONTWRITEBYTECODE=1` set, so no bytecode cache was written for pyro.
- v0.1.0's startup banner prints `torch <environment>` where a version would be expected. This is cosmetic: it comes from the delayed import in `.onLoad`.
- The original v0.2.1 Python environment still exists in the fork checkout (see above). If the v0.1.0 result differs, rerunning in that environment would be a control for the PyPI versus conda builds of torch and numpy. Doing so needs the v0.2.1 R package in some library, and must not write into the checkout.

## Deviations from the brief

- The Python package versions were matched to the original environment instead of being chosen fresh, because that environment was found intact. The brief had said it was gone.
- pytest 9.1.1 was installed because v0.1.0's `_environment.yml` lists it. It is not needed at run time.
- `smoke-test.R` and `smoke-test.log` were added so that the smoke test can be repeated.
- `R CMD INSTALL` was given `--no-multiarch`, which has no effect on macOS arm64.

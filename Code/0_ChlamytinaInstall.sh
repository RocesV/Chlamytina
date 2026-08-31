#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CHLAMYTINA_ENV_NAME:-chlamytina}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$SCRIPT_DIR/.chlamytina-install-logs"
LOG_FILE="$LOG_DIR/install_${STAMP}.log"
WORK_DIR="$(mktemp -d)"

SELFTEST_OUT="$REPO_DIR/tmp/out"
PROT_DP_OUT="$SELFTEST_OUT/proteomics_limma"
PROT_LOLA_OUT="$SELFTEST_OUT/proteomics_limma_lola"
TRANS_DP_OUT="$SELFTEST_OUT/transcriptomics_deseq2_perm"
TRANS_LOLA_OUT="$SELFTEST_OUT/transcriptomics_deseq2_perm_lola"

mkdir -p "$LOG_DIR"
printf "[Chlamytina] Install log started: %s\n" "$STAMP" >"$LOG_FILE"

say() {
  printf "[Chlamytina] %s\n" "$*" | tee -a "$LOG_FILE"
}

say_err() {
  printf "[Chlamytina] %s\n" "$*" | tee -a "$LOG_FILE" >&2
}

cleanup() {
  local rc=$?
  rm -rf "$WORK_DIR"

  if [[ "$rc" -ne 0 ]]; then
    rm -rf "$PROT_DP_OUT" "$PROT_LOLA_OUT" "$TRANS_DP_OUT" "$TRANS_LOLA_OUT"
    rmdir "$SELFTEST_OUT" 2>/dev/null || true
    rmdir "$REPO_DIR/tmp" 2>/dev/null || true
    say_err "Installation failed. Removing environment: $ENV_NAME"
    conda env remove -n "$ENV_NAME" -y >>"$LOG_FILE" 2>&1 || true
    say_err "Command output was saved to the full log, not printed to the terminal."
    say_err "Full log: $LOG_FILE"
  fi

  exit "$rc"
}
trap cleanup EXIT

phase() {
  local label="$1"
  shift
  local start end elapsed rc

  start="$(date +%s)"
  printf "\n" >>"$LOG_FILE"
  say "$label ..."

  if "$@" >>"$LOG_FILE" 2>&1; then
    end="$(date +%s)"
    elapsed=$((end - start))
    say "$label OK (${elapsed}s)"
  else
    rc=$?
    say_err "ERROR in phase: $label"
    say_err "Detailed package output is only in: $LOG_FILE"
    exit "$rc"
  fi
}

env_prefix() {
  conda env list | awk -v env="$ENV_NAME" '$1 == env {print $NF; found=1} END {exit !found}'
}

run_clean_env_bin() {
  local exe="$1"
  shift
  local prefix
  prefix="$(env_prefix)"

  env -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE "$prefix/bin/$exe" "$@"
}

run_clean_env_path() {
  local prefix
  prefix="$(env_prefix)"

  env -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE PATH="$prefix/bin:$PATH" "$@"
}

run_repo_rscript() {
  (
    cd "$REPO_DIR"
    run_clean_env_bin Rscript --vanilla "$@"
  )
}

index_local_channel() {
  local local_channel="$1"

  mkdir -p "$local_channel/noarch" "$local_channel/linux-64"
  rm -rf "$local_channel/.cache" "$local_channel/noarch/.cache" "$local_channel/linux-64/.cache"

  run_clean_env_bin python -m conda_index --threads=1 "$local_channel"
}

check_conda() {
  command -v conda >/dev/null 2>&1
}

check_env_name() {
  case "$ENV_NAME" in
    ""|base|root|*/*|*\\*)
      echo "Unsafe environment name: '$ENV_NAME'"
      return 1
      ;;
  esac
}

remove_existing_env() {
  if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    conda env remove -n "$ENV_NAME" -y
  fi
}

write_local_recipes() {
  local recipes="$WORK_DIR/recipes"
  mkdir -p "$recipes/r-simplecache" "$recipes/r-nvennr" "$recipes/r-cowsay"

  cat > "$recipes/r-simplecache/meta.yaml" <<'YAML'
package:
  name: r-simplecache
  version: 0.5.0

source:
  url:
    - https://cran.r-project.org/src/contrib/simpleCache_0.5.0.tar.gz
    - https://cran.r-project.org/src/contrib/Archive/simpleCache/simpleCache_0.5.0.tar.gz
  sha256: 31fa9e608dcbecf1bb88ca85308d9ada3c41dad911bc297fdf029659550086da

build:
  number: 0
  noarch: generic
  script: $R CMD INSTALL --build .

requirements:
  host:
    - r-base 4.5.3
  run:
    - r-base >=4.5,<4.6.0a0

test:
  commands:
    - $R -e "library('simpleCache')"
YAML

  cat > "$recipes/r-nvennr/meta.yaml" <<'YAML'
package:
  name: r-nvennr
  version: 0.2.3

source:
  url: https://cran.r-project.org/src/contrib/Archive/nVennR/nVennR_0.2.3.tar.gz
  sha256: 6605d2def873552ea372d63cbc2fabaefcebe458312324f006fb0882185daf0c

build:
  number: 0
  script: $R CMD INSTALL --build .
  rpaths:
    - lib/R/lib/
    - lib/

requirements:
  build:
    - {{ compiler('c') }}
    - {{ compiler('cxx') }}
    - make
  host:
    - r-base 4.5.3
    - r-rcpp 1.1.2
  run:
    - r-base >=4.5,<4.6.0a0
    - r-rcpp >=1.1.2

test:
  commands:
    - $R -e "library('nVennR')"
YAML

  cat > "$recipes/r-cowsay/meta.yaml" <<'YAML'
package:
  name: r-cowsay
  version: 1.2.2

source:
  url:
    - https://cran.r-project.org/src/contrib/cowsay_1.2.2.tar.gz
    - https://cran.r-project.org/src/contrib/Archive/cowsay/cowsay_1.2.2.tar.gz
  sha256: 9631bf2085a253bcb7b679bdf155d25f62e0c81e5db55f8e95b354101c2c57b2

build:
  number: 0
  noarch: generic
  script: $R CMD INSTALL --build .

requirements:
  host:
    - r-base 4.5.3
    - r-crayon 1.5.3
    - r-rlang 1.3.0
  run:
    - r-base >=4.5,<4.6.0a0
    - r-crayon >=1.5.3
    - r-rlang >=1.3.0

test:
  commands:
    - $R -e "library('cowsay')"
YAML
}

create_env() {
  conda create -n "$ENV_NAME" -y \
    --no-default-packages \
    --no-pin \
    --override-channels \
    -c conda-forge \
    "python=3.12.*" \
    "mamba=2.9.0" \
    "conda-build=25.7.0" \
    "conda-index=0.12.1"
}

configure_r_library_isolation() {
  local prefix
  prefix="$(env_prefix)"

  mkdir -p "$prefix/etc/conda/activate.d" "$prefix/etc/conda/deactivate.d"

  cat > "$prefix/etc/conda/activate.d/chlamytina-r-libs.sh" <<'SH'
if [ "${R_LIBS+x}" ]; then export CHLAMYTINA_RESTORE_R_LIBS="$R_LIBS"; else unset CHLAMYTINA_RESTORE_R_LIBS; fi
if [ "${R_LIBS_USER+x}" ]; then export CHLAMYTINA_RESTORE_R_LIBS_USER="$R_LIBS_USER"; else unset CHLAMYTINA_RESTORE_R_LIBS_USER; fi
if [ "${R_LIBS_SITE+x}" ]; then export CHLAMYTINA_RESTORE_R_LIBS_SITE="$R_LIBS_SITE"; else unset CHLAMYTINA_RESTORE_R_LIBS_SITE; fi
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE
SH

  cat > "$prefix/etc/conda/deactivate.d/chlamytina-r-libs.sh" <<'SH'
if [ "${CHLAMYTINA_RESTORE_R_LIBS+x}" ]; then export R_LIBS="$CHLAMYTINA_RESTORE_R_LIBS"; else unset R_LIBS; fi
if [ "${CHLAMYTINA_RESTORE_R_LIBS_USER+x}" ]; then export R_LIBS_USER="$CHLAMYTINA_RESTORE_R_LIBS_USER"; else unset R_LIBS_USER; fi
if [ "${CHLAMYTINA_RESTORE_R_LIBS_SITE+x}" ]; then export R_LIBS_SITE="$CHLAMYTINA_RESTORE_R_LIBS_SITE"; else unset R_LIBS_SITE; fi
unset CHLAMYTINA_RESTORE_R_LIBS
unset CHLAMYTINA_RESTORE_R_LIBS_USER
unset CHLAMYTINA_RESTORE_R_LIBS_SITE
SH
}

init_local_channel() {
  local local_channel="$WORK_DIR/local-channel"
  mkdir -p "$local_channel/noarch" "$local_channel/linux-64"
  index_local_channel "$local_channel"
  echo "$local_channel" > "$WORK_DIR/local_channel_path"
}

build_local_package() {
  local recipe_name="$1"
  local recipes="$WORK_DIR/recipes"
  local local_channel
  local_channel="$(cat "$WORK_DIR/local_channel_path")"

  run_clean_env_bin conda-build "$recipes/$recipe_name" \
    --output-folder "$local_channel" \
    --override-channels \
    -c "file://$local_channel" \
    -c conda-forge

  index_local_channel "$local_channel"
}

install_r_stack() {
  local local_channel
  local_channel="$(cat "$WORK_DIR/local_channel_path")"

  run_clean_env_bin mamba install --no-pin -n "$ENV_NAME" -y \
    --override-channels \
    --strict-channel-priority \
    -c "file://$local_channel" \
    -c conda-forge \
    -c bioconda \
    "r-base=4.5.3" \
    "bioconductor-biocversion=3.22.0" \
    "bioconductor-biocgenerics=0.56.0" \
    "bioconductor-biobase=2.70.0" \
    "bioconductor-s4vectors=0.48.0" \
    "bioconductor-iranges=2.44.0" \
    "bioconductor-seqinfo=1.0.0" \
    "bioconductor-genomicranges=1.62.1" \
    "bioconductor-biostrings=2.78.0" \
    "bioconductor-rtracklayer=1.70.1" \
    "bioconductor-annotationdbi=1.72.0" \
    "bioconductor-genomicfeatures=1.62.0" \
    "bioconductor-txdbmaker=1.6.2" \
    "bioconductor-limma=3.66.0" \
    "bioconductor-sva=3.58.0" \
    "bioconductor-lola=1.40.1" \
    "bioconductor-deseq2=1.50.2" \
    "bioconductor-biocparallel=1.44.0" \
    "bioconductor-qvalue=2.42.0" \
    "r-biocmanager=1.30.26" \
    "r-cowsay=1.2.2" \
    "r-optparse=1.8.2" \
    "r-readxl=1.5.0" \
    "r-r.utils=2.13.0" \
    "r-simplecache=0.5.0" \
    "r-nvennr=0.2.3" \
    "r-dplyr=1.2.1" \
    "r-data.table=1.17.8" \
    "r-ggplot2=4.0.3" \
    "r-reshape2=1.4.5" \
    "r-pheatmap=1.0.13" \
    "r-rcolorbrewer=1.1_3" \
    "r-scales=1.4.0" \
    "r-matrixstats=1.5.0" \
    "r-magrittr=2.0.5" \
    "r-rcpp=1.1.2" \
    "r-crayon=1.5.3" \
    "r-rlang=1.3.0"
}

install_jbrowse() {
  run_clean_env_bin mamba install --no-pin -n "$ENV_NAME" -y \
    --override-channels \
    --strict-channel-priority \
    -c conda-forge \
    -c bioconda \
    "jbrowse2=4.3.0"
}

check_r_packages() {
  local check_script="$WORK_DIR/check_r_packages.R"

  cat > "$check_script" <<'RSCRIPT'
if (as.character(getRversion()) != "4.5.3") {
  stop("Unexpected R version: ", as.character(getRversion()), ", expected 4.5.3", call. = FALSE)
}

expected <- c(
  BiocVersion = "3.22.0",
  BiocManager = "1.30.26",
  BiocGenerics = "0.56.0",
  Biobase = "2.70.0",
  S4Vectors = "0.48.0",
  IRanges = "2.44.0",
  Seqinfo = "1.0.0",
  GenomicRanges = "1.62.1",
  Biostrings = "2.78.0",
  rtracklayer = "1.70.1",
  AnnotationDbi = "1.72.0",
  GenomicFeatures = "1.62.0",
  txdbmaker = "1.6.2",
  nVennR = "0.2.3",
  cowsay = "1.2.2",
  optparse = "1.8.2",
  readxl = "1.5.0",
  "R.utils" = "2.13.0",
  simpleCache = "0.5.0",
  limma = "3.66.0",
  sva = "3.58.0",
  LOLA = "1.40.1",
  DESeq2 = "1.50.2",
  BiocParallel = "1.44.0",
  qvalue = "2.42.0",
  dplyr = "1.2.1",
  data.table = "1.17.8",
  ggplot2 = "4.0.3",
  reshape2 = "1.4.5",
  pheatmap = "1.0.13",
  RColorBrewer = "1.1.3",
  scales = "1.4.0",
  matrixStats = "1.5.0",
  magrittr = "2.0.5",
  Rcpp = "1.1.2",
  crayon = "1.5.3",
  rlang = "1.3.0"
)

load_checked <- function(pkg) {
  pkg_path <- system.file(package = pkg)
  if (!nzchar(pkg_path)) {
    stop(
      "R package is not installed in the active R library paths: ", pkg,
      "\n.libPaths():\n", paste(.libPaths(), collapse = "\n"),
      call. = FALSE
    )
  }

  tryCatch(
    suppressPackageStartupMessages(loadNamespace(pkg)),
    error = function(e) {
      stop(
        "R package is installed but failed to load: ", pkg,
        "\nInstalled at: ", pkg_path,
        "\nError: ", conditionMessage(e),
        "\n.libPaths():\n", paste(.libPaths(), collapse = "\n"),
        call. = FALSE
      )
    }
  )

  got <- as.character(utils::packageVersion(pkg))
  if (got != expected[[pkg]]) {
    stop("Unexpected version for ", pkg, ": got ", got, ", expected ", expected[[pkg]], call. = FALSE)
  }

  TRUE
}

invisible(vapply(names(expected), load_checked, logical(1)))

if (!"makeTxDbFromGFF" %in% getNamespaceExports("txdbmaker")) {
  stop("txdbmaker::makeTxDbFromGFF() is not available", call. = FALSE)
}

if (!"transcripts" %in% getNamespaceExports("GenomicFeatures")) {
  stop("GenomicFeatures::transcripts() is not available", call. = FALSE)
}
RSCRIPT

  run_clean_env_bin Rscript "$check_script"
}

check_cli_tools() {
  run_clean_env_path bash -lc '
command -v Rscript >/dev/null
command -v mamba >/dev/null
command -v jbrowse >/dev/null
'
}

check_help() {
  local help_out="$WORK_DIR/1_DataPrepare.help.txt"

  run_clean_env_bin Rscript "$SCRIPT_DIR/1_DataPrepare.R" --help >"$help_out" 2>&1

  grep -q "Usage:" "$help_out"
  grep -q -- "--file1" "$help_out"
  grep -q -- "--condition1" "$help_out"
  grep -q -- "--normalization" "$help_out"
}

require_nonempty_file() {
  local path="$1"

  if [[ ! -s "$path" ]]; then
    echo "Expected non-empty file was not created: $path"
    return 1
  fi
}

require_glob_count() {
  local expected="$1"
  local pattern="$2"
  local files=()

  shopt -s nullglob
  files=( $pattern )
  shopt -u nullglob

  if [[ "${#files[@]}" -ne "$expected" ]]; then
    echo "Expected $expected files matching $pattern, found ${#files[@]}"
    return 1
  fi
}

cleanup_selftest_outputs() {
  rm -rf \
    "$PROT_DP_OUT" \
    "$PROT_LOLA_OUT" \
    "$TRANS_DP_OUT" \
    "$TRANS_LOLA_OUT"

  rmdir "$SELFTEST_OUT" 2>/dev/null || true
  rmdir "$REPO_DIR/tmp" 2>/dev/null || true
}

test_proteomics_dataprepare() {
  run_repo_rscript Code/1_DataPrepare.R \
    -A test/prots/protNitrogen.xlsx -a 4-20 \
    -B test/prots/protCold.xlsx -b 3-15 \
    -C test/prots/protOsmotic.xlsx -q 3-6 \
    -D test/prots/protUV.xlsx -x 3-6 \
    -d TRUE -m limma -n normalizeQuantiles -i TRUE -V v5 \
    -o tmp/out/proteomics_limma
}

check_proteomics_dataprepare_outputs() {
  require_nonempty_file "$PROT_DP_OUT/Diff_background.bed"
  require_nonempty_file "$PROT_DP_OUT/Diff_uniq_file1Treatment2-Treatment1.bed"
  require_nonempty_file "$PROT_DP_OUT/Diff_uniq_file2Treatment2-Treatment1.bed"
  require_nonempty_file "$PROT_DP_OUT/Diff_uniq_file3Treatment2-Treatment1.bed"
  require_nonempty_file "$PROT_DP_OUT/Diff_uniq_file4Treatment2-Treatment1.bed"
}

test_proteomics_lola() {
  run_repo_rscript Code/2_EnrichmentsLOLA.R \
    -A tmp/out/proteomics_limma/Diff_uniq_file1Treatment2-Treatment1.bed \
    -B tmp/out/proteomics_limma/Diff_uniq_file2Treatment2-Treatment1.bed \
    -C tmp/out/proteomics_limma/Diff_uniq_file3Treatment2-Treatment1.bed \
    -D tmp/out/proteomics_limma/Diff_uniq_file4Treatment2-Treatment1.bed \
    -b tmp/out/proteomics_limma/Diff_background.bed \
    -r MMarks -V v5 -m pValueLog -c 4 \
    -o tmp/out/proteomics_limma_lola
}

check_proteomics_lola_outputs() {
  require_nonempty_file "$PROT_LOLA_OUT/Diff_backgroundMMarks_LOLA_results.txt"
  require_nonempty_file "$PROT_LOLA_OUT/Diff_backgroundMMarks.pdf"
}

test_transcriptomics_dataprepare() {
  run_repo_rscript Code/1_DataPrepare.R \
    -A test/trans/Ge_DNMT1mtPlus_WTmtPlus_LengthScaledTPM.txt \
    -a 2-2 \
    -d TRUE -m DESeq2 -n none \
    -f permutations -p 500 -S 123 -c 4 \
    -i TRUE -V v5 -o tmp/out/transcriptomics_deseq2_perm \
    -g TRUE
}

check_transcriptomics_dataprepare_outputs() {
  require_nonempty_file "$TRANS_DP_OUT/file1.bed"
  require_nonempty_file "$TRANS_DP_OUT/file1Treatment2-Treatment1.bed"
  require_nonempty_file "$TRANS_DP_OUT/Diff_background.bed"
  require_nonempty_file "$TRANS_DP_OUT/file1Treatment2-Treatment1.txt"
  require_nonempty_file "$TRANS_DP_OUT/permutations/file1Treatment2-Treatment1/permutation_index.txt"
  require_nonempty_file "$TRANS_DP_OUT/permutations/file1Treatment2-Treatment1/Permutation_001.bed"
  require_nonempty_file "$TRANS_DP_OUT/permutations/file1Treatment2-Treatment1/Permutation_500.bed"
  require_glob_count 500 "$TRANS_DP_OUT/permutations/file1Treatment2-Treatment1/Permutation_*.bed"
}

test_transcriptomics_lola() {
  run_repo_rscript Code/2_EnrichmentsLOLA.R \
    -A tmp/out/transcriptomics_deseq2_perm/file1Treatment2-Treatment1.bed \
    -b tmp/out/transcriptomics_deseq2_perm/file1.bed \
    -p tmp/out/transcriptomics_deseq2_perm/permutations/file1Treatment2-Treatment1 \
    -r CS_Chlamytina -V v5 -m pValueLog -t 1.30103 -e 0.05 -c 4 \
    -o tmp/out/transcriptomics_deseq2_perm_lola
}

check_transcriptomics_lola_outputs() {
  require_nonempty_file "$TRANS_LOLA_OUT/file1CS_Chlamytina_LOLA_results.txt"
  require_nonempty_file "$TRANS_LOLA_OUT/file1CS_Chlamytina_LOLA_permutation_results.txt"
  require_nonempty_file "$TRANS_LOLA_OUT/file1CS_Chlamytina.pdf"
}

phase "Checking conda" check_conda
phase "Checking environment name" check_env_name
phase "Removing existing environment if present" remove_existing_env
phase "Writing local conda recipes" write_local_recipes
phase "Creating $ENV_NAME environment with env-local mamba" create_env
phase "Configuring isolated R library paths for $ENV_NAME" configure_r_library_isolation
phase "Preparing local conda channel with env-local conda-index" init_local_channel

say "Building local R conda packages can take several minutes on the first run: r-simplecache, r-nvennr, r-cowsay."
phase "Building local R conda package r-simplecache=0.5.0" build_local_package r-simplecache
phase "Building local R conda package r-nvennr=0.2.3" build_local_package r-nvennr
phase "Building local R conda package r-cowsay=1.2.2" build_local_package r-cowsay

say "Solving and installing the pinned R 4.5 / Bioconductor 3.22 stack can take a while."
phase "Installing pinned R and Bioconductor stack with mamba" install_r_stack
phase "Installing pinned JBrowse2 with mamba" install_jbrowse
phase "Checking pinned R packages, local R packages, and txdbmaker route" check_r_packages
phase "Checking command-line tools" check_cli_tools
phase "Checking 1_DataPrepare.R help output" check_help

say "Running Chlamytina self-tests with bundled proteomics and transcriptomics data. R output is saved to the install log only."
say "The transcriptomics permutation test runs 500 matched sets and can take several minutes."

phase "Removing previous Chlamytina self-test outputs" cleanup_selftest_outputs
phase "Testing proteomics DataPrepare with limma and normalizeQuantiles" test_proteomics_dataprepare
phase "Checking proteomics DataPrepare BED outputs" check_proteomics_dataprepare_outputs
phase "Testing proteomics LOLA enrichment with MMarks" test_proteomics_lola
phase "Checking proteomics LOLA outputs" check_proteomics_lola_outputs
phase "Testing transcriptomics DataPrepare with DESeq2 permutation background" test_transcriptomics_dataprepare
phase "Checking transcriptomics permutation BED outputs" check_transcriptomics_dataprepare_outputs
phase "Testing transcriptomics LOLA permutation enrichment with CS_Chlamytina" test_transcriptomics_lola
phase "Checking transcriptomics LOLA permutation outputs" check_transcriptomics_lola_outputs
phase "Removing Chlamytina self-test temporary outputs" cleanup_selftest_outputs

say "Installation complete OK"
say "Activate with: conda activate $ENV_NAME"
say "Log saved at: $LOG_FILE"

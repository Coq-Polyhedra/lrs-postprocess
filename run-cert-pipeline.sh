#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="${DATA_DIR:-data}"

LRS_DIR="${LRS_DIR:-/home/labcmap/allamigeon/lrslib-073a}"
LRSGMP="${LRSGMP:-$LRS_DIR/lrsgmp}"

# Build once with:
#   ./run.sh build
#
# Then timed runs use the compiled binary directly, so cargo build/startup
# is not included in the measured time.
BIN="${BIN:-./target/release/lrs-postprocess}"

CERT_GEN="${CERT_GEN:-$BIN postprocess --pretty}"
BIN_CERT_GEN="${BIN_CERT_GEN:-$BIN postprocess --bin}"
CHECKER="${CHECKER:-$BIN check}"

COQC="${COQC:-coqc}"
COQ_TEMPLATE="${COQ_TEMPLATE:-coq/InspectCertificate.v.template}"
COQ_DIR="${COQ_DIR:-${DATA_DIR}/coq}"

now_ms() {
    date +%s%3N
}

format_ms() {
    local ms="$1"

    local s=$((ms / 1000))
    local rem_ms=$((ms % 1000))

    printf "%d.%03d s" "$s" "$rem_ms"
}

run_timed() {
    local name="$1"
    local log_file="$2"
    shift 2

    echo
    echo "=== $name ==="

    local start end elapsed
    start=$(now_ms)

    "$@" > "$log_file" 2>&1

    end=$(now_ms)
    elapsed=$((end - start))

    echo "--- $name completed in $(format_ms "$elapsed") ---"
    echo "log: $log_file"
}

run_timed_stdout_to_file() {
    local name="$1"
    local out_file="$2"
    local log_file="$3"
    shift 3

    echo
    echo "=== $name ==="

    local start end elapsed
    start=$(now_ms)

    "$@" > "$out_file" 2> "$log_file"

    end=$(now_ms)
    elapsed=$((end - start))

    echo "--- $name completed in $(format_ms "$elapsed") ---"
    echo "output: $out_file"
    echo "log: $log_file"
}

build_tools() {
    echo
    echo "=== build Rust tools ==="
    cargo build --release
    echo "--- build completed ---"
}

check_binary() {
    if [[ ! -x "$BIN" ]]; then
        echo "error: compiled binary not found or not executable: $BIN" >&2
        echo "hint: run: $0 build" >&2
        exit 1
    fi
}

compute_ext() {
    local base="$1"
    local ine_file="${DATA_DIR}/${base}.ine"
    local ext_file="${DATA_DIR}/${base}.ext"
    local log_file="${DATA_DIR}/${base}-ext.log"

    if [[ ! -x "$LRSGMP" ]]; then
        echo "error: lrsgmp not found or not executable: $LRSGMP" >&2
        exit 1
    fi

    if [[ ! -f "$ine_file" ]]; then
        echo "error: input .ine file not found: $ine_file" >&2
        exit 1
    fi

    run_timed "compute ext for $base" "$log_file" \
        "$LRSGMP" "$ine_file" "$ext_file"
}

generate_certificate() {
    local base="$1"
    local ine_file="${DATA_DIR}/${base}.ine"
    local ext_file="${DATA_DIR}/${base}.ext"
    local cert_file="${DATA_DIR}/${base}-cert.json"
    local log_file="${DATA_DIR}/${base}-cert.log"

    check_binary

    if [[ ! -f "$ine_file" ]]; then
        echo "error: input .ine file not found: $ine_file" >&2
        exit 1
    fi

    if [[ ! -f "$ext_file" ]]; then
        echo "error: ext file not found: $ext_file" >&2
        exit 1
    fi

    run_timed_stdout_to_file "generate certificate for $base" "$cert_file" "$log_file" \
        $CERT_GEN \
        "$ine_file" \
        "$ext_file"
}

generate_binary_certificate() {
    local base="$1"
    local ine_file="${DATA_DIR}/${base}.ine"
    local ext_file="${DATA_DIR}/${base}.ext"
    local bin_file="${DATA_DIR}/${base}-cert.bin"
    local log_file="${DATA_DIR}/${base}-bin.log"

    check_binary

    if [[ ! -f "$ine_file" ]]; then
        echo "error: input .ine file not found: $ine_file" >&2
        exit 1
    fi

    if [[ ! -f "$ext_file" ]]; then
        echo "error: ext file not found: $ext_file" >&2
        exit 1
    fi

    run_timed_stdout_to_file "generate binary certificate for $base" "$bin_file" "$log_file" \
        $BIN_CERT_GEN \
        "$ine_file" \
        "$ext_file"
}

run_checker() {
    local base="$1"
    local ine_file="${DATA_DIR}/${base}.ine"
    local cert_file="${DATA_DIR}/${base}-cert.json"
    local log_file="${DATA_DIR}/${base}-check.log"

    check_binary

    if [[ ! -f "$ine_file" ]]; then
        echo "error: input .ine file not found: $ine_file" >&2
        exit 1
    fi

    if [[ ! -f "$cert_file" ]]; then
        echo "error: certificate file not found: $cert_file" >&2
        exit 1
    fi

    run_timed "check certificate for $base" "$log_file" \
        $CHECKER \
        "$ine_file" \
        "$cert_file"
}

run_coq_binreader_test() {
    local base="$1"
    local bin_file="${DATA_DIR}/${base}-cert.bin"

    # Coq module names cannot contain '-' characters, so sanitize the generated
    # .v filename. This still keeps the original basename for data/log files.
    local coq_base="${base//-/_}"
    local v_file="${COQ_DIR}/${coq_base}_inspect.v"

    local log_file="${DATA_DIR}/${base}-coq.log"

    if [[ ! -f "$bin_file" ]]; then
        echo "error: binary certificate not found: $bin_file" >&2
        echo "hint: generate it with: $0 bin $base" >&2
        exit 1
    fi

    if [[ ! -f "$COQ_TEMPLATE" ]]; then
        echo "error: Coq template not found: $COQ_TEMPLATE" >&2
        exit 1
    fi

    mkdir -p "$COQ_DIR"

    local abs_bin_file
    abs_bin_file="$(realpath "$bin_file")"

    sed "s|__BIN_FILE__|${abs_bin_file}|g" "$COQ_TEMPLATE" > "$v_file"

    run_timed "Coq/binreader test for $base" "$log_file" \
        "$COQC" "$v_file"
}

clean_generated() {
    local base="$1"

    local ext_file="${DATA_DIR}/${base}.ext"
    local cert_file="${DATA_DIR}/${base}-cert.json"
    local bin_file="${DATA_DIR}/${base}-cert.bin"

    local ext_log="${DATA_DIR}/${base}-ext.log"
    local cert_log="${DATA_DIR}/${base}-cert.log"
    local bin_log="${DATA_DIR}/${base}-bin.log"
    local check_log="${DATA_DIR}/${base}-check.log"
    local coq_log="${DATA_DIR}/${base}-coq.log"

    local coq_base="${base//-/_}"
    local coq_file="${COQ_DIR}/${coq_base}_inspect.v"

    echo
    echo "=== clean generated files for $base ==="

    rm -f \
        "$ext_file" \
        "$cert_file" \
        "$bin_file" \
        "$ext_log" \
        "$cert_log" \
        "$bin_log" \
        "$check_log" \
        "$coq_log" \
        "$coq_file"

    echo "--- cleaned generated files for $base ---"
}

run_one() {
    local cmd="$1"
    local base="$2"

    case "$cmd" in
        ext)
            compute_ext "$base"
            ;;
        cert)
            generate_certificate "$base"
            ;;
        bin)
            generate_binary_certificate "$base"
            ;;
        check)
            run_checker "$base"
            ;;
        coq)
            run_coq_binreader_test "$base"
            ;;
        all)
            compute_ext "$base"
            generate_certificate "$base"
            run_checker "$base"
            generate_binary_certificate "$base"
            run_coq_binreader_test "$base"
            ;;
        clean)
            clean_generated "$base"
            ;;
        *)
            echo "error: unknown command: $cmd" >&2
            exit 1
            ;;
    esac
}

usage() {
    cat <<EOF
Usage:
  $0 build
  $0 ext    BASE [BASE ...]
  $0 cert   BASE [BASE ...]
  $0 bin    BASE [BASE ...]
  $0 check  BASE [BASE ...]
  $0 coq    BASE [BASE ...]
  $0 all    BASE [BASE ...]
  $0 clean  BASE [BASE ...]

For each BASE, the script uses:
  ${DATA_DIR}/BASE.ine
  ${DATA_DIR}/BASE.ext
  ${DATA_DIR}/BASE-cert.json
  ${DATA_DIR}/BASE-cert.bin

Generated logs:
  ${DATA_DIR}/BASE-ext.log
  ${DATA_DIR}/BASE-cert.log
  ${DATA_DIR}/BASE-bin.log
  ${DATA_DIR}/BASE-check.log
  ${DATA_DIR}/BASE-coq.log

Generated Coq files:
  ${COQ_DIR}/BASE_inspect.v
  where '-' in BASE is replaced by '_'.

Examples:
  $0 build
  $0 all cross_9 cross_10 cross_11 cross_12
  $0 ext cross_9
  $0 cert cross_9 cross_10
  $0 bin cross_9
  $0 check cross_12
  $0 coq cross_9
  $0 clean cross_9 cross_10

Environment variables:
  DATA_DIR       directory containing input/output files, default: data
  LRS_DIR        directory containing lrsgmp
  LRSGMP         full path to lrsgmp
  BIN            compiled Rust binary, default: ./target/release/lrs-postprocess
  CERT_GEN       JSON certificate command, default: "\$BIN postprocess --pretty"
  BIN_CERT_GEN   binary certificate command, default: "\$BIN postprocess --bin"
  CHECKER        checker command, default: "\$BIN check"
  COQC           Coq compiler, default: coqc
  COQ_TEMPLATE   Coq template, default: coq/InspectCertificate.v.template
  COQ_DIR        generated Coq files directory, default: data/coq
EOF
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

cmd="$1"
shift

case "$cmd" in
    -h|--help|help)
        usage
        exit 0
        ;;
    build)
        if [[ $# -ne 0 ]]; then
            echo "error: build does not take basenames" >&2
            usage
            exit 1
        fi
        build_tools
        ;;
    ext|cert|bin|check|coq|all|clean)
        if [[ $# -lt 1 ]]; then
            echo "error: expected at least one basename" >&2
            usage
            exit 1
        fi

        for base in "$@"; do
            echo
            echo "########################################"
            echo "# Processing $base"
            echo "########################################"
            run_one "$cmd" "$base"
        done
        ;;
    *)
        echo "error: unknown command: $cmd" >&2
        usage
        exit 1
        ;;
esac
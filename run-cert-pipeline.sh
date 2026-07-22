#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="${DATA_DIR:-data}"

LRS_DIR="${LRS_DIR:-../lrslib-073a}"
LRSGMP="${LRSGMP:-$LRS_DIR/lrsgmp}"
TIME_CMD="${TIME_CMD:-/usr/bin/time}"

# Build once with:
#   ./run-cert-pipeline.sh build
#
# Then timed runs use the compiled binary directly, so cargo build/startup
# is not included in the measured time.
BIN="${BIN:-./target/release/lrs-postprocess}"

JSON_CERT_GEN="${JSON_CERT_GEN:-$BIN postprocess --pretty}"
CERT_GEN="${CERT_GEN:-${BIN_CERT_GEN:-$BIN postprocess --bin}}"

COQC="${COQC:-coqc}"
COQ_TEMPLATE="${COQ_TEMPLATE:-coq/InspectCertificate.v.template}"
COQ_DIR="${COQ_DIR:-${DATA_DIR}/coq}"
LAST_ELAPSED=""

format_elapsed_seconds() {
    local seconds="$1"

    if [[ -z "$seconds" ]]; then
        printf "unknown"
    else
        awk -v t="$seconds" 'BEGIN { printf "%.3f s", t }'
    fi
}

format_elapsed_precise() {
    local seconds="$1"

    if [[ -z "$seconds" ]]; then
        printf "unknown"
    else
        awk -v t="$seconds" 'BEGIN { printf "%.6f s", t }'
    fi
}

format_ratio() {
    local elapsed="$1"
    local reference="$2"

    if [[ -z "$elapsed" || -z "$reference" ]]; then
        printf "unknown"
        return
    fi

    awk -v t="$elapsed" -v r="$reference" 'BEGIN {
        if (r <= 0) printf "unknown";
        else printf "%.3fx lrs", t / r;
    }'
}

format_timing_with_lrs() {
    local elapsed="$1"
    local reference="$2"

    if [[ -z "$elapsed" ]]; then
        printf "unknown"
    elif [[ -z "$reference" ]]; then
        format_elapsed_precise "$elapsed"
    else
        awk -v t="$elapsed" -v r="$reference" 'BEGIN {
            if (r <= 0) printf "%.6f s", t;
            else printf "%.6f s (%.6fx lrs)", t, t / r;
        }'
    fi
}

read_lrs_elapsed() {
    local base="$1"
    local time_file="${DATA_DIR}/${base}-lrs.time"

    if [[ -f "$time_file" ]]; then
        cat "$time_file"
    else
        printf ""
    fi
}

write_lrs_elapsed() {
    local base="$1"
    local elapsed="$2"
    local time_file="${DATA_DIR}/${base}-lrs.time"

    printf "%s\n" "$elapsed" > "$time_file"
}

check_time_command() {
    if [[ ! -x "$TIME_CMD" ]]; then
        echo "error: time command not found or not executable: $TIME_CMD" >&2
        echo "hint: set TIME_CMD to a GNU time-compatible executable" >&2
        exit 1
    fi
}

# Run COMMAND with stdout redirected to OUT_FILE and stderr both displayed live
# and copied to LOG_FILE. Return the command status (or tee's status if the
# command succeeded but logging failed).
run_with_live_stderr() {
    local out_file="$1"
    local log_file="$2"
    shift 2

    local fifo_dir fifo status tee_status tee_pid
    fifo_dir="$(mktemp -d)"
    fifo="${fifo_dir}/stderr"
    mkfifo "$fifo"

    tee "$log_file" < "$fifo" >&2 &
    tee_pid=$!

    if "$@" > "$out_file" 2> "$fifo"; then
        status=0
    else
        status=$?
    fi

    if wait "$tee_pid"; then
        tee_status=0
    else
        tee_status=$?
    fi

    rm -f "$fifo"
    rmdir "$fifo_dir"

    if [[ "$status" -ne 0 ]]; then
        return "$status"
    fi
    return "$tee_status"
}

# Copy the certificate generator's stderr to LOG_FILE and turn its internal
# timing records into the user-facing report on stdout as they arrive. Other
# stderr (warnings and errors) remains on stderr.
stream_certificate_timings() {
    local input_fifo="$1"
    local log_file="$2"
    local reference_seconds="$3"
    local line elapsed label

    : > "$log_file"
    echo "Certificate creation:"

    while IFS= read -r line || [[ -n "$line" ]]; do
        printf '%s\n' "$line" >> "$log_file"

        label=""
        case "$line" in
            "Certificate: read and parse .ine file: "*" s")
                label="read and parse .ine:"
                ;;
            "Certificate: read and parse .ext file: "*" s")
                label="read and parse .ext:"
                ;;
            "Certificate: stream .ext and build vertices/facets: "*" s")
                label="stream .ext and build vertices/facets:"
                ;;
            "Certificate: build vertices and facets: "*" s")
                label="build vertices and facets:"
                ;;
            "Certificate: build facet graph: "*" s")
                label="build facet graph:"
                ;;
            "Certificate: choose root vertex: "*" s")
                label="choose root vertex:"
                ;;
            "Certificate: build root certificate: "*" s")
                label="build root certificate:"
                ;;
            "Certificate: convert inequalities: "*" s")
                label="convert inequalities:"
                ;;
            "Certificate: build geometric graph and lifts: "*" s")
                label="build geometric graph and lifts:"
                ;;
            "Certificate: build full-dimensionality certificate: "*" s")
                label="build full-dimensionality certificate:"
                ;;
            "Certificate: assemble value: "*" s")
                label="assemble certificate value:"
                ;;
            "Create certificate in memory: "*" s")
                label="total certificate creation:"
                ;;
            "Inequality file check: "*" s")
                echo
                echo "Rust certificate check:"
                label="inequality-file check:"
                ;;
            "Well-formedness check: "*" s")
                label="well-formedness check:"
                ;;
            "Uniqueness check: "*" s")
                label="uniqueness check:"
                ;;
            "Feasibility check: "*" s")
                label="feasibility check:"
                ;;
            "Graph check: "*" s")
                label="graph check:"
                ;;
            "Mapping check: "*" s")
                label="mapping check:"
                ;;
            "Root check: "*" s")
                label="root check:"
                ;;
            "Geometric graph check: "*" s")
                label="geometric-graph check:"
                ;;
            "Full dimension check: "*" s")
                label="full-dimension check:"
                ;;
            "Check certificate: "*" s")
                label="total certificate check:"
                ;;
            "Generate and write binary certificate: "*" s")
                echo
                echo "Binary certificate:"
                label="generate and write binary file:"
                ;;
        esac

        if [[ -n "$label" ]]; then
            elapsed="${line% s}"
            elapsed="${elapsed##*: }"
            printf '    %-40s%s\n' "$label" \
                "$(format_timing_with_lrs "$elapsed" "$reference_seconds")"
        elif [[ "$line" != "Generated certificate accepted" ]]; then
            printf '%s\n' "$line" >&2
        fi
    done < "$input_fifo"
}

# Run a certificate command whose stdout is the generated artifact. Internal
# timings are printed on stdout and logged immediately at each phase boundary.
run_with_live_certificate_timings() {
    local out_file="$1"
    local log_file="$2"
    local reference_seconds="$3"
    shift 3

    local fifo_dir fifo status relay_status relay_pid
    fifo_dir="$(mktemp -d)"
    fifo="${fifo_dir}/stderr"
    mkfifo "$fifo"

    stream_certificate_timings "$fifo" "$log_file" "$reference_seconds" &
    relay_pid=$!

    if "$@" > "$out_file" 2> "$fifo"; then
        status=0
    else
        status=$?
    fi

    if wait "$relay_pid"; then
        relay_status=0
    else
        relay_status=$?
    fi

    rm -f "$fifo"
    rmdir "$fifo_dir"

    if [[ "$status" -ne 0 ]]; then
        return "$status"
    fi
    return "$relay_status"
}

# Timed command whose stdout/stderr both go to the log file.
# Arguments:
#   run_timed NAME LOG_FILE REFERENCE_SECONDS COMMAND...
# If REFERENCE_SECONDS is nonempty, the report includes the ratio to lrs time.
run_timed() {
    local name="$1"
    local log_file="$2"
    local reference_seconds="$3"
    shift 3

    check_time_command

    echo
    echo "=== $name ==="

    local time_file status elapsed ratio
    time_file="$(mktemp)"

    if "$TIME_CMD" -f "%e" -o "$time_file" "$@" > "$log_file" 2>&1; then
        status=0
    else
        status=$?
    fi

    elapsed="$(cat "$time_file" 2>/dev/null || true)"
    rm -f "$time_file"

    ratio="$(format_ratio "$elapsed" "$reference_seconds")"

    if [[ "$status" -eq 0 ]]; then
        if [[ -n "$reference_seconds" ]]; then
            echo "--- $name completed in $(format_elapsed_seconds "$elapsed") (${ratio}) ---"
        else
            echo "--- $name completed in $(format_elapsed_seconds "$elapsed") ---"
        fi
    else
        if [[ -n "$reference_seconds" ]]; then
            echo "--- $name failed after $(format_elapsed_seconds "$elapsed") (${ratio}) ---"
        else
            echo "--- $name failed after $(format_elapsed_seconds "$elapsed") ---"
        fi
    fi
    echo "log: $log_file"

    LAST_ELAPSED="$elapsed"
    return "$status"
}

# Timed command whose stdout is an output file and stderr is a log file.
# Arguments:
#   run_timed_stdout_to_file NAME OUT_FILE LOG_FILE REFERENCE_SECONDS COMMAND...
# If REFERENCE_SECONDS is nonempty, the report includes the ratio to lrs time.
run_timed_stdout_to_file() {
    local name="$1"
    local out_file="$2"
    local log_file="$3"
    local reference_seconds="$4"
    shift 4

    check_time_command

    echo
    echo "=== $name ==="

    local time_file status elapsed ratio
    time_file="$(mktemp)"

    if run_with_live_stderr "$out_file" "$log_file" \
        "$TIME_CMD" -f "%e" -o "$time_file" "$@"; then
        status=0
    else
        status=$?
    fi

    elapsed="$(cat "$time_file" 2>/dev/null || true)"
    rm -f "$time_file"

    ratio="$(format_ratio "$elapsed" "$reference_seconds")"

    if [[ "$status" -eq 0 ]]; then
        if [[ -n "$reference_seconds" ]]; then
            echo "--- $name completed in $(format_elapsed_seconds "$elapsed") (${ratio}) ---"
        else
            echo "--- $name completed in $(format_elapsed_seconds "$elapsed") ---"
        fi
    else
        if [[ -n "$reference_seconds" ]]; then
            echo "--- $name failed after $(format_elapsed_seconds "$elapsed") (${ratio}) ---"
        else
            echo "--- $name failed after $(format_elapsed_seconds "$elapsed") ---"
        fi
    fi
    echo "output: $out_file"
    echo "log: $log_file"

    LAST_ELAPSED="$elapsed"
    return "$status"
}

# Run a certificate command whose stdout is the generated artifact.
run_stdout_to_file() {
    local name="$1"
    local out_file="$2"
    local log_file="$3"
    local reference_seconds="$4"
    shift 4

    echo
    echo "=== $name ==="

    local status
    if run_with_live_certificate_timings \
        "$out_file" "$log_file" "$reference_seconds" "$@"; then
        status=0
    else
        status=$?
    fi

    if [[ "$status" -eq 0 ]]; then
        echo "--- $name completed ---"
    else
        echo "--- $name failed ---"
    fi
    echo "output: $out_file"
    echo "log: $log_file"

    return "$status"
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

    # lrs is the reference step, so its ratio is 1.000x lrs.
    run_timed "compute ext for $base" "$log_file" "" \
        "$LRSGMP" "$ine_file" "$ext_file"

    write_lrs_elapsed "$base" "$LAST_ELAPSED"
    echo "reference lrs time: $(format_elapsed_seconds "$LAST_ELAPSED") (1.000x lrs)"
}

generate_certificate() {
    local base="$1"
    local ine_file="${DATA_DIR}/${base}.ine"
    local ext_file="${DATA_DIR}/${base}.ext"
    local cert_file="${DATA_DIR}/${base}-cert.json"
    local log_file="${DATA_DIR}/${base}-cert.log"
    local lrs_elapsed
    lrs_elapsed="$(read_lrs_elapsed "$base")"

    check_binary

    if [[ ! -f "$ine_file" ]]; then
        echo "error: input .ine file not found: $ine_file" >&2
        exit 1
    fi

    if [[ ! -f "$ext_file" ]]; then
        echo "error: ext file not found: $ext_file" >&2
        exit 1
    fi

    if [[ -z "$lrs_elapsed" ]]; then
        echo "warning: missing lrs time for $base; run '$0 lrs $base' first for ratios" >&2
    fi

    run_timed_stdout_to_file "generate certificate for $base" "$cert_file" "$log_file" "$lrs_elapsed" \
        $JSON_CERT_GEN \
        "$ine_file" \
        "$ext_file"
}

generate_binary_certificate() {
    local base="$1"
    local ine_file="${DATA_DIR}/${base}.ine"
    local ext_file="${DATA_DIR}/${base}.ext"
    local bin_file="${DATA_DIR}/${base}-cert.bin"
    local log_file="${DATA_DIR}/${base}-bin.log"
    local lrs_elapsed
    lrs_elapsed="$(read_lrs_elapsed "$base")"

    check_binary

    if [[ ! -f "$ine_file" ]]; then
        echo "error: input .ine file not found: $ine_file" >&2
        exit 1
    fi

    if [[ ! -f "$ext_file" ]]; then
        echo "error: ext file not found: $ext_file" >&2
        exit 1
    fi

    if [[ -z "$lrs_elapsed" ]]; then
        echo "warning: missing lrs time for $base; run '$0 lrs $base' first for ratios" >&2
    fi

    run_stdout_to_file "create, check, and encode certificate for $base" "$bin_file" "$log_file" "$lrs_elapsed" \
        $CERT_GEN \
        "$ine_file" \
        "$ext_file"
}

run_coq_binreader_test() {
    local base="$1"
    local bin_file="${DATA_DIR}/${base}-cert.bin"

    # Coq module names cannot contain '-' characters, so sanitize the generated
    # .v filename. This still keeps the original basename for data/log files.
    local coq_base="${base//-/_}"
    local v_file="${COQ_DIR}/${coq_base}_inspect.v"

    local log_file="${DATA_DIR}/${base}-coq.log"
    local lrs_elapsed
    lrs_elapsed="$(read_lrs_elapsed "$base")"

    if [[ ! -f "$bin_file" ]]; then
        echo "error: binary certificate not found: $bin_file" >&2
        echo "hint: generate it with: $0 cert $base" >&2
        exit 1
    fi

    if [[ ! -f "$COQ_TEMPLATE" ]]; then
        echo "error: Coq template not found: $COQ_TEMPLATE" >&2
        exit 1
    fi

    if [[ -z "$lrs_elapsed" ]]; then
        echo "warning: missing lrs time for $base; run '$0 lrs $base' first for ratios" >&2
    fi

    mkdir -p "$COQ_DIR"

    local abs_bin_file
    abs_bin_file="$(realpath "$bin_file")"

    sed "s|__BIN_FILE__|${abs_bin_file}|g" "$COQ_TEMPLATE" > "$v_file"

    run_timed "Coq/binreader test for $base" "$log_file" "$lrs_elapsed" \
        "$COQC" "$v_file"
}

clean_generated() {
    local base="$1"

    local ext_file="${DATA_DIR}/${base}.ext"
    local cert_file="${DATA_DIR}/${base}-cert.json"
    local bin_file="${DATA_DIR}/${base}-cert.bin"
    local lrs_time_file="${DATA_DIR}/${base}-lrs.time"

    local ext_log="${DATA_DIR}/${base}-ext.log"
    local cert_log="${DATA_DIR}/${base}-cert.log"
    local bin_log="${DATA_DIR}/${base}-bin.log"
    local coq_log="${DATA_DIR}/${base}-coq.log"

    local coq_base="${base//-/_}"
    local coq_file="${COQ_DIR}/${coq_base}_inspect.v"

    echo
    echo "=== clean generated files for $base ==="

    rm -f \
        "$ext_file" \
        "$cert_file" \
        "$bin_file" \
        "$lrs_time_file" \
        "$ext_log" \
        "$cert_log" \
        "$bin_log" \
        "$coq_log" \
        "$coq_file"

    echo "--- cleaned generated files for $base ---"
}

run_one() {
    local cmd="$1"
    local base="$2"

    case "$cmd" in
        lrs|ext)
            compute_ext "$base"
            ;;
        json)
            generate_certificate "$base"
            ;;
        cert|bin)
            generate_binary_certificate "$base"
            ;;
        rocq|coq)
            run_coq_binreader_test "$base"
            ;;
        run|all)
            compute_ext "$base"
            generate_binary_certificate "$base"
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
  $0 PATTERN[.ine] [PATTERN[.ine] ...]
  $0 run    PATTERN[.ine] [PATTERN[.ine] ...]
  $0 build
  $0 lrs    PATTERN[.ine] [PATTERN[.ine] ...]
  $0 cert   PATTERN[.ine] [PATTERN[.ine] ...]
  $0 rocq   PATTERN[.ine] [PATTERN[.ine] ...]
  $0 clean  PATTERN[.ine] [PATTERN[.ine] ...]

The default command is 'run'. It executes the complete production pipeline:
  lrs -> create certificate -> Rust check -> binary encoding -> Rocq check

Individual production stages:
  lrs     Generate BASE.ext with lrsgmp.
  cert    From an existing BASE.ext, create and Rust-check BASE-cert.bin.
  rocq    Check an existing BASE-cert.bin with Rocq.

JSON diagnostics (not needed by the production pipeline):
  $0 json   PATTERN[.ine] [PATTERN[.ine] ...]

Compatibility aliases:
  all = run, ext = lrs, bin = cert, coq = rocq

Instance selection:
  Every PATTERN is a Bash extended regular expression matched against the
  complete basename of each ${DATA_DIR}/*.ine file, without the .ine extension.
  Thus 'cross10' selects exactly cross10.ine, while 'cross(8|9|10)' selects
  cross8.ine, cross9.ine, and cross10.ine. A trailing .ine is accepted.
  Matches are processed in lexicographic order and duplicates are removed.
  The command fails if a pattern is invalid or matches no file.

For each BASE, the script uses:
  ${DATA_DIR}/BASE.ine
  ${DATA_DIR}/BASE.ext
  ${DATA_DIR}/BASE-cert.json
  ${DATA_DIR}/BASE-cert.bin
  ${DATA_DIR}/BASE-lrs.time

Generated logs:
  ${DATA_DIR}/BASE-ext.log
  ${DATA_DIR}/BASE-cert.log
  ${DATA_DIR}/BASE-bin.log
  ${DATA_DIR}/BASE-coq.log

Generated Coq files:
  ${COQ_DIR}/BASE_inspect.v
  where '-' in BASE is replaced by '_'.

Examples:
  $0 build
  $0 cross_9 cross_10 cross_11 cross_12
  $0 lrs cross_9.ine
  $0 cert cross_9 cross_10
  $0 cert 'cross(8|9|10)'
  $0 'dual_cyclic_d(1[4-9]|20)_n.*'
  $0 rocq cross_9
  $0 json cross_12
  $0 clean cross_9 cross_10

Timing:
  The ext step, i.e. lrsgmp, is used as the reference lrs time.
  The binary pipeline reports every certificate-construction phase, every
  Rust checker group, and binary generation/writing as soon as it completes.
  Each value is printed immediately in seconds and as a multiple of the
  measured lrs time; the raw timing records are also retained in BASE-bin.log.
  The reference is stored in ${DATA_DIR}/BASE-lrs.time.
  The run command uses the production path and does not generate JSON:
  lrs -> in-memory Rust check -> binary encoding -> Rocq check.

Environment variables:
  DATA_DIR       directory containing input/output files, default: data
  LRS_DIR        directory containing lrsgmp
  LRSGMP         full path to lrsgmp
  TIME_CMD       GNU time-compatible command, default: /usr/bin/time
  BIN            compiled Rust binary, default: ./target/release/lrs-postprocess
  CERT_GEN       binary certificate command, default: "\$BIN postprocess --bin"
  JSON_CERT_GEN  diagnostic JSON command, default: "\$BIN postprocess --pretty"
  BIN_CERT_GEN   deprecated fallback name for CERT_GEN
  COQC           Coq compiler, default: coqc
  COQ_TEMPLATE   Coq template, default: coq/InspectCertificate.v.template
  COQ_DIR        generated Coq files directory, default: data/coq
EOF
}

expand_patterns() {
    local pattern_arg pattern anchored path filename base regex_status
    local -a all_bases=() matches=()
    local -A selected=()

    if [[ ! -d "$DATA_DIR" ]]; then
        echo "error: data directory not found: $DATA_DIR" >&2
        exit 1
    fi

    while IFS= read -r path; do
        filename="${path##*/}"
        all_bases+=("${filename%.ine}")
    done < <(find "$DATA_DIR" -maxdepth 1 -type f -name '*.ine' -print | LC_ALL=C sort)

    for pattern_arg in "$@"; do
        pattern="${pattern_arg%.ine}"
        anchored="^(${pattern})$"
        matches=()

        if [[ "" =~ $anchored ]]; then
            :
        else
            regex_status=$?
            if [[ "$regex_status" -eq 2 ]]; then
                echo "error: invalid regular expression: $pattern_arg" >&2
                exit 1
            fi
        fi

        for base in "${all_bases[@]}"; do
            if [[ "$base" =~ $anchored ]]; then
                matches+=("$base")
                selected["$base"]=1
            fi
        done

        if [[ ${#matches[@]} -eq 0 ]]; then
            echo "error: pattern matched no .ine basenames in $DATA_DIR: $pattern_arg" >&2
            exit 1
        fi
    done

    for base in "${all_bases[@]}"; do
        if [[ -n "${selected[$base]+x}" ]]; then
            printf '%s\n' "$base"
        fi
    done
}

run_patterns() {
    local command="$1"
    shift
    local expanded base
    local -a bases

    expanded="$(expand_patterns "$@")"
    mapfile -t bases <<< "$expanded"

    for base in "${bases[@]}"; do
        echo
        echo "########################################"
        echo "# Processing $base"
        echo "########################################"
        run_one "$command" "$base"
    done
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
    run|all|lrs|ext|cert|bin|rocq|coq|json|clean)
        if [[ $# -lt 1 ]]; then
            echo "error: expected at least one instance pattern" >&2
            usage
            exit 1
        fi
        run_patterns "$cmd" "$@"
        ;;
    *)
        # With no explicit command, treat every argument as an instance pattern
        # and run the complete production pipeline on all matching basenames.
        set -- "$cmd" "$@"
        run_patterns run "$@"
        ;;
esac

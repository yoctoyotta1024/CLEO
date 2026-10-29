#! /bin/bash

zarr_name=$1 # e.g. "fromfile_sol.zarr" or "fromfile_irreg_sol.zarr"
PYTHON=${PYTHON:-python} # python with numpy and zarr, used to report size of differences

zarr1=./build/bin_1/${zarr_name}
if [ ! -d "$zarr1" ]; then
  echo "$zarr1 does not exist."
  exit 1
fi

# print the largest difference between the data of each named array in two zarr datasets
report_max_diffs() {
  "$PYTHON" - "$@" <<'EOF'
import sys
import numpy as np
import zarr

zarr1, zarrX = sys.argv[1], sys.argv[2]
for name in sys.argv[3:]:
    try:
        a = np.asarray(zarr.open(f"{zarr1}/{name}", mode="r")[...])
        b = np.asarray(zarr.open(f"{zarrX}/{name}", mode="r")[...])
    except Exception as e:
        print(f"    {name}: could not read array ({e})")
        continue

    if a.shape != b.shape:
        print(f"    {name}: shapes differ, {a.shape} vs {b.shape}")
        continue

    if not (np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number)):
        ndiff = np.count_nonzero(a != b)
        print(f"    {name}: {ndiff} of {a.size} elements differ (non-numeric dtype)")
        continue

    a64, b64 = a.astype(np.float64), b.astype(np.float64)
    nan_a, nan_b = np.isnan(a64), np.isnan(b64)
    differs = ~((a64 == b64) | (nan_a & nan_b))
    ndiff = np.count_nonzero(differs)
    if ndiff == 0:
        print(f"    {name}: data identical (only metadata or encoding differs)")
        continue

    absdiff = np.where(nan_a | nan_b, 0.0, np.abs(a64 - b64))
    idx = np.unravel_index(np.argmax(absdiff), absdiff.shape)
    maxdiff = absdiff[idx]
    scale = max(abs(a64[idx]), abs(b64[idx]))
    reldiff = maxdiff / scale if scale > 0 else float("inf")
    nnan = np.count_nonzero(nan_a != nan_b)
    print(f"    {name}: {ndiff} of {a.size} elements differ, "
          f"max abs diff = {maxdiff:.6g} (rel {reldiff:.3g}) at index {tuple(int(i) for i in idx)}: "
          f"{a[idx]} vs {b[idx]}" + (f", {nnan} NaN mismatches" if nnan else ""))
EOF
}

status=0
for processes in 2 4 8; do
    zarrX=./build/bin_${processes}/${zarr_name}
    if [ ! -d "$zarrX" ]; then
        echo "$zarrX does not exist."
        exit 1
    fi

    echo "--- Comparing sequential run with ${processes} processes ---"
    failed=()
    run_failed=0

    # group-level metadata
    for meta in .zattrs .zgroup; do
        if diff -q "$zarr1/$meta" "$zarrX/$meta" > /dev/null 2>&1; then
            echo "  PASS  $meta"
        else
            echo "  FAIL  $meta"
            diff "$zarr1/$meta" "$zarrX/$meta" | sed 's/^/        /'
            run_failed=1
        fi
    done

    # each array: all chunks plus its .zarray and .zattrs
    for array in $(ls "$zarr1/"); do
        if diff -rq "$zarr1/$array" "$zarrX/$array" > /dev/null 2>&1; then
            echo "  PASS  $array"
        else
            echo "  FAIL  $array"
            diff -rq "$zarr1/$array" "$zarrX/$array" 2>&1 | sed 's/^/        /'
            failed+=("$array")
        fi
    done

    if [ ${#failed[@]} -gt 0 ]; then
        run_failed=1
        echo "  Size of differences in failed arrays:"
        report_max_diffs "$zarr1" "$zarrX" "${failed[@]}"
    fi

    if [ $run_failed -ne 0 ]; then
      status=1
      echo "Run with ${processes} processes has different results than sequential run"
    fi
done;

if [ $status -ne 0 ]; then
  exit 1
fi

echo "All parallel execution results match the sequential run results"
exit 0;

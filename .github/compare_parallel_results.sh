#! /bin/bash

for processes in 2 4 8; do
    if [ -s "diffs" ]; then
        rm diffs
    fi

    for file in $(ls bin_1/fromfile_sol.zarr/); do
        filename=$(ls bin_1/fromfile_sol.zarr/$file)
        hexdump bin_1/fromfile_sol.zarr/$file/$filename > original
        hexdump bin_${processes}/fromfile_sol.zarr/$file/$filename > new
        diff original new >> diffs
    done;
    rm new original

    diff bin_1/fromfile_sol.zarr/.zattrs bin_${processes}/fromfile_sol.zarr/.zattrs >> diffs
    diff bin_1/fromfile_sol.zarr/.zgroup bin_${processes}/fromfile_sol.zarr/.zgroup >> diffs
    for file in $(ls bin_1/fromfile_sol.zarr/); do
        filename=$(ls -la bin_1/fromfile_sol.zarr/$file)
        diff bin_1/fromfile_sol.zarr/$file/.zattrs bin_${processes}/fromfile_sol.zarr/$file/.zattrs >> diffs
        diff bin_1/fromfile_sol.zarr/$file/.zarray bin_${processes}/fromfile_sol.zarr/$file/.zarray >> diffs
    done;

    if [ -s "diffs" ]; then
      echo "Run with ${processes} processes has different results than sequential run"
      exit 1;
    fi
done;

echo "All parallel execution results match the sequential run results"
exit 0;

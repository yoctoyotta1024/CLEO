#! /bin/bash

zarr_name=$1 # e.g. "fromfile_sol.zarr" or "fromfile_irreg_sol.zarr"

zarr1=./build/bin_1/${zarr_name}
if [ ! -d $zarr1 ]; then
  echo "$zarr1 does not exist."
  exit 1
fi

for processes in 2 4 8; do
    if [ -s "diffs" ]; then
        rm diffs
    fi

    zarrX=./build/bin_${processes}/${zarr_name}
    if [ ! -d $zarrX ]; then
        echo "$zarrX does not exist."
        exit 1
    fi

    for file in $(ls $zarr1/); do
        filename=$(ls $zarr1/$file)
        hexdump $zarr1/$file/$filename > original
        hexdump $zarrX/$file/$filename > new
        diff original new >> diffs
    done;
    rm new original

    diff $zarr1/.zattrs $zarrX/.zattrs >> diffs
    diff $zarr1/.zgroup $zarrX/.zgroup >> diffs
    for file in $(ls $zarr1/); do
        filename=$(ls -la $zarr1/$file)
        diff $zarr1/$file/.zattrs $zarrX/$file/.zattrs >> diffs
        diff $zarr1/$file/.zarray $zarrX/$file/.zarray >> diffs
    done;

    if [ -s "diffs" ]; then
      echo "Run with ${processes} processes has different results than sequential run"
      exit 1;
    fi
done;

echo "All parallel execution results match the sequential run results"
exit 0;

# Pairwise clustering by SNV sets
This script takes a table of sample IDs and VCF file paths and outputs an SNV difference matrix and cluster assignments according to the SNV difference threshold defined by the user.

## Inputs
You can view options by running `Rscript pairwiseclustering.r --help`.
```
usage: pairwiseclustering.r [-h] -i IDS -f FRONTCLIP -b BACKCLIP [-o ODIR]
                            [-t THRESHOLD]

optional arguments:
  -h, --help                          show this help message and exit
  -i IDS, --ids IDS                   file with sample IDs
  -f FRONTCLIP, --frontclip FRONTCLIP bases to clip from the front
  -b BACKCLIP, --backclip BACKCLIP    bases to clip from the back
  -o ODIR, --odir ODIR                output directory
  -t THRESHOLD, --threshold THRESHOLD SNV difference threshold for clustering
```

| Flag | Description                                                                                                                  |
| ---- | ---------------------------------------------------------------------------------------------------------------------------- |
| `-i` | CSV file with two columns. The left column should contain sample IDs and the right column should contain paths to VCF files. |
| `-o` | Output directory. The default is `$PWD`.                                                                                     |
| `-f` | The number of nucleotides to clip from the front end of the reference sequence.                                              |
| `-b` | The number of nucleotides to clip from the back end of the reference sequence.                                               |
| `-t` | The SNV error threshold to use when assigning clusters. The default is `0` SNVs.                                             |

## Outputs
The script outputs the difference matrix and cluster assignments into the output directory of your choice.

| Filename          | Contents                                                                  |
| ----------------- | ------------------------------------------------------------------------- |
| `differences.csv` | Difference matrix. Rows and columns are in the same order.                |
| `clusters.csv`    | Table of samples IDs (left column) and numeric cluster IDs (right column) |

# Kmer_spectrum_Hw3
## Description
K-mer spectrum and genome size estimator in Python.
## Theoretical principle
This tool relies on k-mer based genome size estimation. That is, k-mer spectrum  is parsed in 'number-frequency' coordinates; then, the highest peak (except for noise region clinging to ordinate axis) is used as a reference and each number tag is multiplied by its frequency; finally, sum of such multiplications is divided by chosen *k* (length of k-mer) value.
## Code description
The code body comprises of two main parts:
- *Class **kmer_spectrun***: a class object standing for kmer spectrum for input file;
- an *argparse* CL integrative object.

### *Class **kmer_spectrum***
This class possesses several functions that perform rough genome size estimation:
- **read_walking** and **library_walking**: two functions that enable kmer registration over all reads in an inpu file;
- **array_builder**: builds a *pandas* array comprising all found number-frequency pairs;
- **noise_eraser**: an **extremely** rough noise cutter that finds minimum in the first half of frequency distribution and uses it as a cutting baseline for future genome size estimation;
- **genome_size_est**: gives a rough assessment of genome size following aformentioned logic (baseline positive).
- **visualize**: draws a kmer number-frequency distribution.
A **finalize** function creates a **kmer_spectrum** object and launches all needed functions one after another.

### *argparse*
The CL arguments used are:
- **-i**: an input *fastq* file;
- **-k**: a size of kmer;
- **-q**: a quality baseline; all kmrs with PHRED quality beyond this minimum will be omitted;
- **-o**: a name for output *png* picture

## Example

Estimated genome size in this case is

## Acknowledgements
Eugene Bakin of Bioinformatics Institute for his Python crash course

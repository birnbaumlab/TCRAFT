# TCRAFT
TCRAFT is a a scalable and cost-effective platform for pooled assembly of synthetic TCR libraries. This software package converts a table of user-inputted TCR sequences into CDR3 oligonucleotide pools to order for the TCRAFT assembly process. 

The TCRAFT package requires prior installation of Python 3.7 or higher as well as an updated version of pip.

## Installation
### Option 1: Quick Install from GitHub
1. Open Terminal or Command Prompt on your machine.
2. Enter the following command: `pip install git+https://github.com/birnbaumlab/TCRAFT.git`

### Option 2: Build the TCRAFT package locally
1. Download or clone the TCRAFT git repository locally.
2. In Terminal or Command Prompt, navigate to the base directory of the package (`/TCRAFT`) and execute the following command: `pip install .`
  
### Backup: Install Dependencies Manually
By default, pip will automatically find and install all Python packages that are necessary for TCRAFT to run. In the event that dependency errors arise during installation, you can also manually install the requisite packages in your Python environment by running the following line in Terminal or Command Prompt:

`pip install numpy pandas tqdm biopython`

TCRAFT uses the following package versions:
- numpy>=1.26.4
- pandas>=2.2.3
- tqdm>=4.67.0
- biopython>=1.84

## Running TCRAFT
Once installed, TCRAFT can be executed from the terminal using two simple commands: `TCRAFT-generate` and `TCRAFT-validate`. `TCRAFT-generate` takes in a CSV-formatted list of TCR sequences to assemble and generates lists of oligonucleotide sequences to order in pooled format for the TCRAFT assembly process. `TCRAFT-validate` takes in the list of oligonucleotides created by `TCRAFT-generate` and simulates TCRAFT Golden Gate assembly of these oligos into complete TCR sequences. `TCRAFT-validate` then checks all of these simulated assemblies to ensure that the desired TCR sequence is correctly assembled and does not contain any errors.

### Example Pipeline:

1. Prepare a CSV table of TCR sequences for TCRAFT assembly. The CSV table **must** contain the following columns or else `TCRAFT-generate` will throw an error:
    - `V_alpha` (TCR $\alpha$ chain V region IMGT ID)
    - `V_beta`  (TCR $\beta$ chain V region IMGT ID)
    - `J_alpha` (TCR $\alpha$ chain J region IMGT ID)
    - `J_beta` (TCR $\beta$ chain J region IMGT ID)
    - `CDR3_alpha` (TCR $\alpha$ chain CDR3 amino acid sequence)
    - `CDR3_beta` (TCR $\beta$ chain CDR3 amino acid sequence)

    A sample TCR list is provided in this repository as an example of how to properly format and prepare the TCR list and can be used to test the TCRAFT pipeline.

2. Generate CDR3 oligo pools by executing the following command: `TCRAFT-generate [name of TCR list file].csv --output_dir [output directory path]`

    The `--output_dir` flag is optional and allows you to specify the output file directory. If not specified the output files will be stored in `./TCRAFT-generate_<current date>`. 

    `TCRAFT-generate` creates the following output files in the output directory:
    - `All_<Number of TCRs>_CDR3_oligos.csv`: Contains all successfully generated oligos
    - `Under300_All_CDR3_oligos.csv`: Contains only oligos with length <= 300bp (which can be ordered from suppliers such as Twist Biosciences).
    - `Over300_All_CDR3_oligos.csv`: Contains only oligos with length > 300bp (which can be ordered from alternative suppliers such as IDT).
    - `Under300_Pool<A or B>_CDR3_oligos.csv`: Split up the under 300bp oligos into pools that can be ordered together from an oligo pool supplier.
    - `Over300_Pool<A or B>_CDR3_oligos.csv`: Split up the oversize >300bp oligos into pools that can be ordered together from an oligo pool supplier.
    - `Oligo_Pool_Metadata.csv`: Table of oligo pool statistics which are helpful to plan out TCRAFT experimental protocols (ex: number of oligos per pool, average oligo size, etc.). This table can be directly pasted into the TCRAFT experiment planning workbook template to automatically compute resuspension volumes, mixing ratios, etc.
    - `Invalid_TCR_Table.csv`: Contains TCRs that failed TCRAFT oligo generation. The most common reason for failure is due to specification of a V/J allele ID which is not supported by TCRAFT (TCRAFT only uses functional, coding alleles as annotated by IMGT).

3. Validate the CDR3 oligo pools by executing the following command: `TCRAFT-validate [name of TCRAFT oligo pool file].csv --output_dir [output directory path]`

    Similar to above, the `--output_dir` flag is optional. If not specified the output files will be stored in `./TCRAFT-validate_<current date>`. 

    `TCRAFT-validate` will simulate Golden Gate assembly of all inputted CDR3 oligos. If an oligo fails assembly validation, either due to presence of extraneous restriction sites or mismatch between the simulated and desired TCR protein sequence, the program will print the oligo/TCR that needs fixing. At the end, the program will output the percentage of inputted TCRs that passed assembly validation.

    `TCRAFT-validate` creates the following output files in the output directory:
    - `Reference_Assembled_TCR_Sequences.csv`: These are the full successfully validated TCR sequences obtained after `TCRAFT-validate` simulated TCRAFT assembly of the CDR3 oligos into full TCR sequences. This reference list is sometimes useful to cross-check downstream experimental sequencing data against the ground-truth list of assembled sequences.
    - `Failed_Assemblies.csv`: Contains the TCRs that failed assembly validation along with info that is helpful for debugging the failures.
    - `Assembly_Metadata.csv`: Similar to the `Oligo_Pool_Metadata.csv` file, this table contains useful statistics for experimental planning and can also be directly pasted into the TCRAFT experiment planning workbook template.



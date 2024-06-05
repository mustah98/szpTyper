# szpTyper

**`szpTyper`** is a Python tool designed for genotypic typing of *Streptococcus equi* genomes based on variants from the *SzP* (Streptococcus zooepidemicus M-like protein) gene. Leveraging NCBI's BLAST and Databases and BioPython libraries, **`szpTyper`** automates the process of querying and analyzing *SzP* gene sequences against a provided *S. equi* genome, aiding in the identification of specific *SzP* variants. Furhtermore, this tool alwoes for automated downlaod and updating of all certain *SzP* vartiants currently available on NCBIs database. 

## Database
szpTyper currently inlcudes **201** distinct *SzP* genes. 

## Dependencies

- Python 3.x
- BioPython
- NCBI BLAST

## Installation

You can install **`szpTyper`** using `git`:

```bash
git clone https://github.com/your_username/szpTyper.git
cd szpTyper
````

## Usage **`szp_DB_download`**

```bash
szp_DB_download.py 
```
This will fetch all SzP Sequences that fulfill a sequence length of 900-1500 bp, download them in order to create a indexed database for querying. 

## Usage **`szpTyper`** Single File

```bash
szpTyper.py -i <input_genome_file.fasta> -o <output_directory> -f <szp_genes_fasta> -t <num_threads>
```

## Usage **`szpTyper`** Multiple Files
```bash 
szpTyper.py -i <file directory> -o <output_directory> --summary -f <szp_genes_fasta> -t <num_threads>
```

With input variables being: 

    -i, --input: Path to the input genome file or directory.
    -o, --output: Output directory for the results.
    -f, --fasta: Path to the SzP genes FASTA file.
    -t, --threads: Number of threads for BLASTN (default: 4).
    -s, --summary: Creates a summary of all results if multiple genomes are   scaned.



## Further Descriptions
SzP Gene: The *SzP* gene, equivalent to the Streptococcus zooepimdicus M-Like protein gene, plays a crucial role in streptococcal infections and pathogenesis. szpTyper aids in understanding the genetic diversity of S. equi strains by identifying and characterizing SzP gene variants.

BLASTN Analysis: szpTyper utilizes NCBI's BLASTN tool to compare *SzP* gene sequences against the provided genome, determining the closest matches and their characteristics.

Molecular Epidemiology: With its ability to rapidly genotype S. equi strains, szpTyper facilitates molecular epidemiology studies, enabling researchers to track outbreaks, investigate transmission dynamics, and monitor genetic variations over time.

Veterinary Research: szpTyper is invaluable in veterinary research, aiding in the surveillance and control of S. equi infections in horses and other animals. Its high-throughput capability accelerates strain characterization and informs vaccination strategies.

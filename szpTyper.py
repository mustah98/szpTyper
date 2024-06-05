from Bio import SeqIO
import argparse
import os
import subprocess

def check_database(db_name):
    """ Checks if the BLAST database exists and is correctly specified in the input path. """
    db_prefix = os.path.join(db_name, "SzP_sequences")  # Adjusting to the specific naming convention of database files
    if os.path.exists(f"{db_prefix}.nin"):
        return True
    else:
        print(f"BLAST database '{db_name}' does not exist or is not correctly specified. Expected database files with prefix 'SzP_sequences'.")
        return False

def get_gene_lengths(fasta_file):
    """ Returns a dictionary of gene lengths keyed by gene IDs. """
    gene_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract the genomic szp_type from the sequence ID
        szp_type = record.id.split('_')[0]
        gene_lengths[record.id] = (len(record.seq), szp_type)
    return gene_lengths

def run_blast(query_file, db_name, output_file, num_threads, gene_lengths):
    """ Runs BLASTN, calculates coverage, and appends genomic szp_types. """
    db_prefix = os.path.join(db_name, "SzP_sequences")  # Ensure we are using the correct prefix for the database
    headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    temp_output = output_file + ".tmp"
    cmd = (
        f"blastn -query {query_file} -db {db_prefix} -out {temp_output} "
        f"-outfmt '6 {headers}' -max_target_seqs 1 -num_threads {num_threads}"
    )
    subprocess.run(cmd, shell=True, check=True)

    # Calculate coverage and append results, and find the entry with the highest coverage
    highest_coverage_entry = None
    highest_coverage_value = 0
    with open(temp_output, "r") as f_temp, open(output_file, "w") as f_out:
        f_out.write(headers.replace(' ', '\t') + '\tcoverage\tSzP_Type\n')  # Append headers for coverage and genomic szp_type
        for line in f_temp:
            parts = line.strip().split()
            sseqid = parts[1]
            sstart = int(parts[8])
            send = int(parts[9])
            gene_length, szp_type = gene_lengths[sseqid]
            coverage = abs(send - sstart + 1) / gene_length * 100
            summary_line = line.strip() + f'\t{coverage:.2f}\t{szp_type}'
            f_out.write(summary_line + '\n')
            # Track the highest coverage
            if coverage > highest_coverage_value:
                highest_coverage_value = coverage
                highest_coverage_entry = [os.path.basename(query_file), parts[0], sseqid, szp_type, parts[2], f'{coverage:.2f}']
    os.remove(temp_output)
    return highest_coverage_entry

def main():
    parser = argparse.ArgumentParser(description="SzP Typing Tool")
    parser.add_argument("-i", "--input", help="Input genome file or directory")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-d", "--db", help="Path to the BLAST database with 'SzP_Database' prefix")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("--summary", action="store_true", help="Generate a summary file of all results")
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if not args.db.startswith("SzP_Database"):
        print("Database path must start with 'SzP_Database'.")
        return

    if not check_database(args.db):
        return

    if os.path.isfile(args.input):
        genome_files = [args.input]
    elif os.path.isdir(args.input):
        genome_files = [os.path.join(args.input, file) for file in os.listdir(args.input) if file.endswith((".fna", ".fasta"))]
    else:
        print("Invalid input. Please provide a valid genome file or directory.")
        return
    
    gene_lengths = get_gene_lengths(os.path.join(args.db, "SzP_sequences.fasta"))

    all_summary_data = []
    for genome_file in genome_files:
        output_file = os.path.join(args.output, os.path.splitext(os.path.basename(genome_file))[0] + "_results.tsv")
        highest_coverage_entry = run_blast(genome_file, args.db, output_file, args.threads, gene_lengths)
        if highest_coverage_entry:
            all_summary_data.append(highest_coverage_entry)
        
        print(f"Results saved to {output_file}")

    if args.summary:
        summary_output_file = os.path.join(args.output, "summary_results.tsv")
        with open(summary_output_file, "w") as f:
            f.write("Filename\tqseqid\tsseqid\tSzP_Type\tpident\tcoverage\n")
            for data in all_summary_data:
                f.write(f"{data[0]}\t{data[1]}\t{data[2]}\t{data[3]}\t{data[4]}\t{data[5]}\n")
        print(f"Summary results saved to {summary_output_file}")

if __name__ == "__main__":
    main()

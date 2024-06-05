from Bio import Entrez, SeqIO
from datetime import datetime
import os
import subprocess

def parse_date(record):
    """ Attempt to parse the date from a GenBank record. """
    try:
        # Try to extract and parse the date from the annotations
        return datetime.strptime(record.annotations['date'], "%d-%b-%Y")
    except KeyError:
        # If the date is missing or the field has unexpected formatting, use a fallback date
        return datetime.min  # This places records with unknown dates at the beginning

def create_blast_db(fasta_file, db_name, output_dir):
    """ Create a BLAST database from a FASTA file. """
    cmd = f"makeblastdb -in {fasta_file} -dbtype nucl -out {os.path.join(output_dir, db_name)}"
    subprocess.run(cmd, shell=True, check=True)
    print("BLAST database created successfully.")

def fetch_sequences(email, output_dir, min_length=900):
    """ Fetch sequences from NCBI, sort by submission date, add an index, save to a FASTA file, and create a BLAST database. """
    Entrez.email = email  # Always provide your email address
    fasta_filename = "SzP_sequences.fasta"
    output_fasta = os.path.join(output_dir, fasta_filename)
    db_name = os.path.splitext(fasta_filename)[0]

    # Define the query for the SzP gene in Streptococcus equi
    query = '(Streptococcus SzP AND 900:1500[Sequence Length])'

    # Use Entrez.esearch to search for gene IDs
    search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=10000)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results['IdList']

    if not ids:
        print("No sequences found.")
        return

    # Fetch the sequences with detailed record information
    fetch_handle = Entrez.efetch(db="nucleotide", id=','.join(ids), rettype="gb", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "gb"))
    fetch_handle.close()

    # Sort records by parsed submission date in ascending order
    records.sort(key=parse_date, reverse=False)

    # Add index to each record's description
    for index, record in enumerate(records, start=1):
        record.id = f"{index}_{record.id}"
        record.description = f"{index} {record.description}"

    # Write all records to a file
    SeqIO.write(records, output_fasta, "fasta")
    print(f"{len(records)} sequences written to {output_fasta}")

    # Create BLAST database
    create_blast_db(output_fasta, db_name, output_dir)

if __name__ == "__main__":
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"./SzP_Database_{current_time}"  # Adjust as needed
    user_email = "email@example.com"  # Replace with your actual email

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fetch_sequences(user_email, output_dir)

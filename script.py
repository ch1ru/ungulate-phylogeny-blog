import shutil
from typing import List
import pandas as pd
import os
import json
import subprocess
import numpy as np
from tqdm import tqdm
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
from Bio.Nexus import Nexus

def download_fasta(genbank_file: str, out_dir: str, format = "fasta"):

    # make folder
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    with open(genbank_file) as f:
        data = json.load(f)

    os.chdir(out_dir) # set output directory for subprocess

    length = len(data['accession'])
    print(f"Collecting {length} items from Genbank...")

    for item in tqdm(data['accession']):
        if item.find(':') >= 0:
            first, last = item.split(':')[0], item.split(':')[1]
            # strip letter prefix
            prefix = first[:2]
            first_number, last_number = int(first[2:]), int(last[2:])
            for i in tqdm(np.arange(start=first_number, stop=last_number, step=1), leave=False):
                item = prefix + str(i)
                if os.path.isfile(f"{item}.fa"):
                    continue
                    
                subprocess.call(["ncbi-acc-download", "--format", format, item.lstrip()])

        if os.path.isfile(f"{item}.fa"):
            continue
            
        subprocess.call(["ncbi-acc-download", "--format", format, item.lstrip()])
    
    os.chdir('..')


def concat_seqs(in_dir: str, out_dir: str):
    seq_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    else:
        shutil.rmtree(out_dir) # since append mode is in use, we don't want to duplicate records, so delete previous ones
        os.makedirs(out_dir)

    for file in seq_files:
        for seq_record in SeqIO.parse(os.path.join(in_dir, file), "fasta"):
            record: SeqRecord = seq_record
            species = f"{record.description.split(' ')[1]}_{record.description.split(' ')[2]}"
            if species == "Erinaceus_europaeus":
                species = "Erinaceus"
            reg = re.search(r"\([^\(\)]*\)\sgene.", record.description)
            exon = re.search(r"exon [1-9]\d*", record.description)
            if not reg:
                continue
            gene = reg.group().lstrip().removeprefix('(').split(')')[0]
            if exon:
                gene = f"{gene}_{exon.group()}"

            # remove -s and spaces from name and make upper case
            gene = gene.replace(' ', '_').replace('-','').upper()

            # write to new fasta file, reassign id as species name, contains all records for that gene

            with open(os.path.join(out_dir, f"{gene}.fasta"), "a") as f: # append mode
                new_record = record
                new_record.id = species
                f.write(f">{new_record.id} {new_record.description}\n")
                f.write(f"{str(new_record.seq)}\n")


def align_sequences(in_dir, out_dir, outfmt="nexus"):

    # make folder
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    aligned_files = []
    unaligned_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]

    for file in unaligned_files:
        outfile = f"{file.split('.')[0]}.{'fa' if outfmt == 'fasta' else 'nex'}"
        cmd = ["clustalw", "-type=dna", f"-infile={in_dir + '/' + file}", f"-output={outfmt}", f"-outfile={out_dir}/{outfile}", "-outorder=input"]
        if os.path.exists(outfile):
            if input(f"{outfile} already exists, overwrite? [y/N] ") == "y": # if exists, ask user if they want to replace it
                subprocess.run(cmd)
        
        else: # no alignments exist, proceed with multiple alignment
            subprocess.run(cmd)
        
        aligned_files.append(outfile)

    return aligned_files


def concat_nexus(in_dir, output_name: str):
    aligned_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]

    os.chdir(in_dir)
  
    nexi = [(fname, Nexus.Nexus(fname)) for fname in aligned_files]

    combined = Nexus.combine(nexi)
    allowed_exts = ["nex", "nxs", "nexus"]
    output_name = f"{output_name}.nex" if not any(x in output_name for x in allowed_exts) else output_name

    os.chdir('..')
    combined.write_nexus_data(filename=open(output_name, "w"))


def integrate(in_dir: str, into: str, records: List[str]):

    seq_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]
    original = [f for f in listdir(into) if isfile(join(into, f))]

    for file in seq_files:
        for seq_record in SeqIO.parse(os.path.join(in_dir, file), "fasta"):
            if file in original:
                if seq_record.id in records:
                    with open(os.path.join(into, file), "a") as f: # append mode
                        f.write(f">{seq_record.id} {seq_record.description}\n")
                        f.write(f"{str(seq_record.seq.replace('-', ''))}\n")


def copy_files(source_dir: str, target_dir: str):
    file_names = os.listdir(source_dir)
    
    for file_name in file_names:
        shutil.copy(os.path.join(source_dir, file_name), target_dir)



if __name__ == "__main__":

    download_fasta(genbank_file="genbank.json", out_dir='tmp_seqs')

    concat_seqs(in_dir='tmp_seqs', out_dir='seqs')

    copy_files(source_dir='extras', target_dir='seqs')

    integrate(in_dir='orm2', into='seqs', records=['Equus', 'Bos', 'Erinaceus', 'Sorex', 'Vicugna', 'Canis', 'Felis' ])

    align_sequences(in_dir='seqs', out_dir='aligned_seqs')

    concat_nexus(in_dir='aligned_seqs', output_name="partition")

    
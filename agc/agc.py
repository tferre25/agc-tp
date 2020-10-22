#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
#import nwalign3 as nw

__author__ = "Théo Ferreira"
__copyright__ = "Universite de Paris"
__credits__ = ["Théo Ferreira"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Théo Ferreira"
__email__ = "theo.ferreira.med@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

# ===========================================================================
#                   1.  Dé-duplication en séquence “complète”
# ===========================================================================

def read_fasta(amplicon_file, minseqlen):
    '''
    Function to read a .fasta.gz file
    args :
        amplicon_file: str
        minseqlen: int
    Returns a generator of sequences of size >=  minseqlen
    '''
    if amplicon_file.endswith(".gz"):
        fillin = gzip.open(amplicon_file, "rb")
    else:
        fillin = open(amplicon_file)

    fasta_seq = ""
    for line in fillin:
        if line.startswith('>'):

            if len(fasta_seq) >= minseqlen:
                yield fasta_seq
            fasta_seq = ""
            
        else :
            fasta_seq += line.strip()
    yield fasta_seq
    fillin.close()

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    '''
    Function to
    args :
    Returns :
    '''
    sequence_list = [seq for seq in read_fasta(amplicon_file, minseqlen)]
    
    # Return a list of the n most common elements and their counts from the most common to the least.
    # If n is omitted or None, most_common() returns all elements in the counter.
    # Elements with equal counts are ordered arbitrarily.
    for sequence in Counter(sequence_list).most_common():
        if sequence[1] >= mincount:
            yield sequence

def get_chunks(sequence, chunk_size):
    '''
    Function to
    args :
    Returns :
    '''
    chunk_seq_list = []

    for i in range(0, len(sequence), chunk_size):
        chunk = sequence[i:i + chunk_size]
        if len(chunk) == chunk_size:
            chunk_seq_list.append(sequence[i:i + chunk_size])
    
        if len(chunk_seq_list) >= 4:
            return chunk_seq_list

def cut_kmer(sequence, kmer_size):
    '''
    Function to cut the sequence in kmer with a k size
    args :
        sequence : sequences in fastq
        kmer_size : size of kmer
    returns a generator of kmers of the sequence
    '''
    for letter in range(len(sequence) - kmer_size + 1):
        yield sequence[letter: letter + kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    '''
    Function to
    args :
    Returns :
    '''
    for seq in cut_kmer(sequence, kmer_size):
        if seq not in kmer_dict:
            kmer_dict[seq] = [id_seq]
        elif id_seq not in kmer_dict[seq]:
            kmer_dict[seq].append(id_seq)
    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    '''
    Function to
    args :
    Returns :
    '''
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]

def get_identity(alignment_list):
    '''
    Function to
    args :
    Returns :
    '''
    count_same = 0

    for i in range(len(alignment_list[0])) :
        if alignment_list[0][i] == alignment_list[1][i] :
            count_same =+ 1
    return count_same/len(alignment_list[0]) * 100

def detect_chimera(perc_identity_matrix):
    '''
    Merci Yann pour les explications et pour les lignes de code
    '''
    list_std = []
    booleen_seq0 = False
    booleen_seq1 = False

    for id_list in perc_identity_matrix : 
        list_std.append(statistics.stdev(id_list))
        if id_list[0] > id_list[1] : 
            booleen_seq0 = True
        if id_list[0] < id_list[1] :  
            booleen_seq1 = True

    if statistics.mean(list_std) > 5 and booleen_seq1 and booleen_seq0:
        return True
    else :
        return False

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    '''
    Function to
    args :
    Returns :
    '''

    # kmer_dict = {}
    # no_chim_list = []

    # dereplication = dereplication_fulllength(amplicon_file, minseqlen, mincount)


    # for sequence in dereplication : 
    #     seq_list = get_chunks(sequence, kmer_size)


    #     for seq in seq_list :
    #         chunck = search_mates(kmer_dict,seq, kmer_size)

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    '''
    Function to
    args :
    Returns :
    '''    
    # OTU_list
    # chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)

def fill(text : str, width = 80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file : str):
    '''
    A écrire en binaire ...
    '''    
    with open(output_file, "w") as fillout:
        for i, OTU in enumerate(OTU_list):
            fillout.write(f">OTU_{i + 1}, occurence : {OTU[1]}" + "\n")
            fillout.write(fill(OTU[0]))
            fillout.write("\n")

# ==============================================================
# Main program
# ==============================================================


def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()

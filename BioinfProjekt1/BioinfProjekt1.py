import sys
import subprocess
import os

import Bio
from Bio.Alphabet import generic_dna
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.SubsMat import FreqTable

alph=Gapped(IUPAC.ambiguous_dna)

def printAlignmentInfo(alignment, alphabet):
    seqlist = []
    for record in alignment:
        seqlist.append(record.seq)               

    m = motifs.create(seqlist, alphabet)
    pwm = m.counts.normalize()
    consensus = pwm.consensus

    print(alignment)

    print('first description: %s' % alignment[0].description)
    print('first sequence: %s' % alignment[0].seq)
    print('matrix %s' % pwm)
    print('consensus %s' % consensus)    
    return

def loadAlignmentFromFile( fileName, format ):
    alignment = AlignIO.read(fileName, format, alphabet = alph)    
    return alignment

def loadAlignmentsListFromFile(fileName, format):
    alignments = list(AlignIO.parse(fileName, format))
    return alignments

def alignSequences(filename):
    clustalw_exe = r"clustalw2.exe"
    cline = ClustalwCommandline(clustalw_exe, infile=filename, outfile='test.aln', gapopen = 0, gapext = 0)

    return_code = subprocess.call(str(cline), shell=(sys.platform != "win32"))
    assert return_code == 0, "Calling ClustalW failed"

    resultAlignment = loadAlignmentFromFile('test.aln', 'clustal')
    return resultAlignment

def combineAlignments(alignments):
    AlignIO.write(alignments, "temp.fasta", "fasta")

    result = alignSequences("temp.fasta")    
    return result

option = input("Single alignment (s) or multiple (m)? ")
if(option == "s" or option == "S"):
    fileName = input("Enter filename ")
    format = input("Enter format ")
    alignment = loadAlignmentFromFile(fileName, format)
    printAlignmentInfo(alignment, alph)
elif(option == "m" or option == "M"):
    fileName = input("Enter filename ")
    format = input("Enter format ")
    alignments = loadAlignmentsListFromFile(fileName, format)
    result = combineAlignments(alignments)
    printAlignmentInfo(result, alph)
elif(option == "d" or option == "D"):
    alignments = loadAlignmentsListFromFile("my_example.phy", "phylip")
    result = combineAlignments(alignments)
    printAlignmentInfo(result, alph)
else:
    print("Unrecognized option")

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

    summary_align = AlignInfo.SummaryInfo(alignment)

    consensus2 = summary_align.dumb_consensus()
    my_pssm = summary_align.pos_specific_score_matrix(consensus,
                                                  chars_to_ignore=['N'])

    print(alignment)

    print('first description: %s' % alignment[0].description)
    print('first sequence: %s' % alignment[0].seq)
    print('length %i' % alignment.get_alignment_length())


    print('matrix pwm %s' % pwm)
    print('consensus (motifs) %s' % consensus)
    
    print('matrix pssm %s' % my_pssm)
    print('consensus (AlignInfo.SummaryInfo) %s' % consensus2) 

    return

def loadAlignmentFromFile( fileName, format ):
    alignment = AlignIO.read(fileName, format, alphabet = alph)    
    return alignment

def loadAlignmentsListFromFile(fileName, format):
    alignments = list(AlignIO.parse(fileName, format))
    return alignments

def performAlignSequences(filename):
    clustalw_exe = r"clustalw2.exe"
    cline = ClustalwCommandline(clustalw_exe, infile=filename, outfile='alignOutput.aln', gapopen = 0, gapext = 0)
    return_code = subprocess.call(str(cline), shell=(sys.platform != "win32"))
    assert return_code == 0, "Calling ClustalW failed"

    resultAlignment = loadAlignmentFromFile('alignOutput.aln', 'clustal')
    return resultAlignment

def alignProfiles(profile1, profile2):
    clustalw_exe = r"clustalw2.exe"
    cline = ClustalwCommandline(clustalw_exe, profile1=profile1, profile2=profile2, outfile='alignOutput.aln', gapopen = 0, gapext = 0, profile = True)
    return_code = subprocess.call(str(cline), shell=(sys.platform != "win32"))
    assert return_code == 0, "Calling ClustalW failed"

    resultAlignment = loadAlignmentFromFile('alignOutput.aln', 'clustal')
    return resultAlignment


def alignSequences(alignments):
    AlignIO.write(alignments, "temp.fasta", "fasta")

    result = performAlignSequences("temp.fasta")    
    return result

print("Wybierz opcję:")
print("1. wyświetlenie informacji o wielodopasowaniu")
print("2. złożenie dwóch wielodopasowań")
print("3. progressive multialigning")

option = input()
if(option == "1"):
    fileName = input("Podaj nazwę pliku z wielodopasowaniem ")
    format = input("Podaj format ")
    alignment = loadAlignmentFromFile(fileName, format)
    printAlignmentInfo(alignment, alph)
elif(option == "2"):
    profile1 = input("Podaj nazwę pliku z 1. wielodopasowaniem ")
    profile2 = input("Podaj nazwę pliku z 2. wielodopasowaniem ")
    result = alignProfiles(profile1, profile2)
    printAlignmentInfo(result, alph)
elif(option == "3"):
    fileName = input("Podaj nazwę pliku z sekwencjami ")
    format = input("Podaj format ")
    alignments = loadAlignmentsListFromFile(fileName, format)
    result = alignSequences(alignments)
    printAlignmentInfo(result, alph)
elif(option == "1d"):
    alignment = loadAlignmentFromFile("prof1.fasta", "fasta")
    printAlignmentInfo(alignment, alph)
elif(option == "2d"):
    result = alignProfiles("prof1.fasta", "prof2.fasta")
    printAlignmentInfo(result, alph)
elif(option == "3d"):
    alignments = loadAlignmentsListFromFile("my_example.phy", "phylip")
    result = alignSequences(alignments)
    printAlignmentInfo(result, alph)
else:
    print("Unrecognized option")

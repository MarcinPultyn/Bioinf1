import Bio
from Bio.Alphabet import generic_dna
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC, Gapped

#instances = [Seq("GGCGTTCAGGCA"),
#             Seq("AAGAATCAGTCA"),
#             Seq("CAAGGAGTTCGC"),
#             Seq("CACGTCAATCAC"),
#             Seq("CAATAATATTCG"),
#            ]

#m = motifs.create(instances)
#pwm = m.counts.normalize()

#print(pwm)
alph=Gapped(IUPAC.unambiguous_dna)
align1 = MultipleSeqAlignment([
             SeqRecord(Seq("ACTGCTAGCTAG", alph), id="Alpha"),
             SeqRecord(Seq("ACT-CTAGCTAG", alph), id="Beta"),
             SeqRecord(Seq("ACTGCTAGCTAG", alph), id="Gamma"),
         ], alph)

align2 = MultipleSeqAlignment([
             SeqRecord(Seq("GTCAGC-AG", Gapped(generic_dna)), id="Delta"),
             SeqRecord(Seq("GACAGCTAG", Gapped(generic_dna)), id="Epsilon"),
             SeqRecord(Seq("GTCAGCTAG", Gapped(generic_dna)), id="Zeta"),
         ], Gapped(generic_dna))

align3 = MultipleSeqAlignment([
             SeqRecord(Seq("ACTAGTACAGCTG", Gapped(generic_dna)), id="Eta"),
             SeqRecord(Seq("ACTAGTACAGCT-", Gapped(generic_dna)), id="Theta"),
             SeqRecord(Seq("-CTACTACAGGTG", Gapped(generic_dna)), id="Iota"),
         ], Gapped(generic_dna))

#align4 = MultipleSeqAlignment([
#             SeqRecord(Seq("ABC-A", generic_dna), id="Eta"),
#             SeqRecord(Seq("ABABA", generic_dna), id="Theta"),
#             SeqRecord(Seq("ACCB-", generic_dna), id="Iota"),
#             SeqRecord(Seq("CB_BC", generic_dna), id="Zeta"),
#         ])

my_alignments = [align1, align2, align3]

seqlist = []
for record in align1:
    seqlist.append(record.seq)

m = motifs.create(seqlist, alph)
pwm = m.counts.normalize()
print(pwm)
consensus = pwm.consensus
print(consensus)

#summary_align = AlignInfo.SummaryInfo(align4)
#consensus = summary_align.dumb_consensus()

#my_pssm = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore = ['X'])

#print(consensus)
#print(my_pssm)


AlignIO.write(my_alignments, "my_example.fasta", "fasta")
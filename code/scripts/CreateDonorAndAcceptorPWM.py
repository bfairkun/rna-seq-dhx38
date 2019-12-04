#!/usr/bin/env python3

from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
import sys
import Bio

DonorsFile, AcceptorsFile, DonorsToScore, AcceptorsToScore, DonorsOut, AcceptorsOut = sys.argv[1:]

print("Reading annotated splice donors")
Donors = []
with open(DonorsFile, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        Donors.append(record.seq.upper())
DonorMotif = motifs.create(Donors)

print("Creating position specific scoring matrix")
donor_pssm = DonorMotif.counts.normalize().log_odds()

print("Scoring provided donors and writing to {}".format(DonorsOut))
DonorsOut_fh=open(DonorsOut, "w")
with open(DonorsToScore, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        mySeq = Seq(str(record.seq.upper()), Bio.Alphabet.IUPAC.unambiguous_dna)
        DonorsOut_fh.write('{0}\t{1}\t{2:.4g}\n'.format(record.id.split("(")[0], mySeq, donor_pssm.calculate(mySeq)))
DonorsOut_fh.close()




print("Reading annotated splice acceptors")

Acceptors = []
with open(AcceptorsFile, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        Acceptors.append(record.seq.upper())
AcceptorMotif = motifs.create(Acceptors)

print("Creating position specific scoring matrix")
acceptor_pssm = AcceptorMotif.counts.normalize().log_odds()

print("Scoring provided acceptors and writing to {}".format(AcceptorsOut))
AcceptorsOut_fh=open(AcceptorsOut, "w")
with open(AcceptorsToScore, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        mySeq = Seq(str(record.seq.upper()), Bio.Alphabet.IUPAC.unambiguous_dna)
        AcceptorsOut_fh.write('{0}\t{1}\t{2:.4g}\n'.format(record.id.split("(")[0], mySeq, acceptor_pssm.calculate(mySeq)))
AcceptorsOut_fh.close()


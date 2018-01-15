"""
-----
Description : Align a set of sequence pairwise calulate their similarity and return it in a Matrix
-----
Author : Marie BLUNTZER
-----
Version 0.2
-----
"""
import matplotlib.pyplot as plt
import numpy as np
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
from Bio import pairwise2
from matplotlib import cm as cm
import argparse


def pairwiseblossum(sequences):
    matrix = matlist.blosum62
    scores = np.zeros((len(sequences),len(sequences)),dtype=float)
    for i in range(len(sequences)):
        for j in range(len(sequences)):
         A=pairwise2.align.globaldx(sequences[i].seq, sequences[j].seq, matrix)
         if i==j :
             sc=A[0][2]
         scores[i,j]=A[0][2]/sc*100
    return scores


def pairwise(sequences):
    scores = np.zeros((len(sequences),len(sequences)),dtype=float)
    for i in range(len(sequences)):
        for j in range(len(sequences)):
         A= pairwise2.align.localmx(sequences[i].seq, sequences[j].seq, 1 , 0)
         scores[i,j]=A[0][2]/min(len(sequences[i].seq), len(sequences[j].seq)) *100
    return scores


def displayresults(scores,labels):
    cmap = cm.get_cmap('Reds')
    fig, ax = plt.subplots(figsize=(15,15))
    cax = ax.matshow(scores, interpolation='nearest', cmap=cmap)
    ax.grid(True)
    plt.xticks(range(scores.shape[0]), labels, rotation=90);
    plt.yticks(range(scores.shape[0]), labels);
    fig.colorbar(cax,)

    for i in range(scores.shape[0]):
        for j in range(scores.shape[0]):
          ax.annotate(int(scores[i,j]),xy=(i-0.25,j))
    plt.show()

def main():

    parser = argparse.ArgumentParser(description='Produce a similarity matrix')
    parser.add_argument('-ff', '--fastafile' , required=True, type=str, help='fasta file containing the sequences to compare')
    parser.add_argument('-m', default='binary', choices=['binary','blossum'] , help='method to compare the sequences : if binary option the score will be the pourcentage of similar residus. if blossum option the score will be calculated using a blossum matrix : similar residu will score more than totaly different ones')
    args = parser.parse_args()
    print(vars(args))
    sequences=tuple(SeqIO.parse(args.fastafile, "fasta"))
    print ("The following %s sequences will be compared pairwise " %len(sequences) )
    labels = []
    for sequence in SeqIO.parse(args.fastafile ,"fasta") :
        print(sequence.id)
        labels.append(sequence.id.split("|")[-1])
    if args.m=='binary':
        scores=pairwise(sequences)
    elif args.m=='blossum':
        scores=pairwiseblossum(sequences)
    displayresults(scores,labels)

if __name__ == "__main__":
    main()

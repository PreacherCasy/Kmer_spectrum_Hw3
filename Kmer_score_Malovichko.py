import argparse
import collections
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO


class Kmer_spectrum:
    def __init__(self, file):
        self.file = file
        self.k = 1
        self.q = 1
        self.gb = collections.Counter()
        self.grid = pd.DataFrame({'Number':[0], 'Frequency':[0]})
        self.G = int()
        self.baseline = int()

    def library_walking(self, k, q):
        self.k, self.q = k, q

        def read_walking(read, k, q):
            seq_lng = len(read)
            kmer_list = []
            for index in range(seq_lng - k + 1):
                if not (kmer_list is None):
                    counter = 0
                    for value in read.letter_annotations['phred_quality'][index:index + k]:
                        if value >= q:
                            counter += 1
                        else:
                            break
                    if counter == k:
                        kmer = str(read.seq[index:index + k])
                        kmer_list.append(kmer)
                else:
                    kmer = str(read.seq[index:index + k])
                    kmer_list.append(kmer)
            return kmer_list

        self.global_count = collections.defaultdict(int)
        for read in SeqIO.parse(self.file, 'fastq'):
            Kmer_spector = read_walking(read, self.k, self.q)
            for kmer in Kmer_spector:
                self.global_count[kmer] += 1

    def array_builder(self):
        bg = collections.Counter(sorted(list(self.global_count.values())))
        self.grid = pd.DataFrame({'Number': sorted([j for j in bg.keys()]), 'Frequency':[i for i in bg.values()]})

    def noise_eraser(self):
        noise_region = list(self.grid['Frequency'])[0:len(self.grid.index)//2]
        self.baseline = min(noise_region)*2

    def genome_size_est(self, baseline=False):
        if baseline:
            self.baseline = baseline
        else:
            self.noise_eraser()
        numerator = sum([list(self.grid['Number'])[self.baseline:][i] * list(self.grid['Frequency'])[self.baseline:][i] for i in range(len(self.grid.index[self.baseline:]))])
        self.G = round(numerator / self.k)

    def visualize(self, ax='auto'):

        sns.set(style="white", palette="BuGn_r")

        plt.title(f'Kmer spectrum (k={self.k}, q={self.q})')

        df = self.grid
        df['Frequency'].plot(logy=True)
        plt.ylabel('Frequency')
        plt.xlabel('Number of kmers')

        plt.axvline(x=self.baseline, color='r')

        plt.axis(ax)
        if o:
            plt.savefig(o + '.png')

        plt.show()

def finalize(file, k, q, b):
    spectrum = Kmer_spectrum(file)
    spectrum.library_walking(k, q)
    spectrum.array_builder()
    spectrum.genome_size_est(b)
    spectrum.visualize()
    return spectrum.G, '\n', spectrum.baseline



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A fast k-mer score counter with built-in approximate genome size estimator')
    parser.add_argument('-i', '--input', help='Path to input fastq file', metavar='Str', type=str)
    parser.add_argument('-k', '--kmer', help='A length of analized k-mers', metavar = 'Int', type = int, default=15)
    parser.add_argument('-q', '--quality', help='Per base quality baseline for k-mers. Sequences with at least any base with quality below this value will ber dropped out from'
                                   'the analysis', metavar = 'Int', type = int, default=35)
    parser.add_argument('-b', '--baseline', help='Custom cutoff baseline. If none, autofunction result will be used', metavar='Int', type=int)
    parser.add_argument('-o', '--output' , help='Output fastq file', nargs='?', metavar='Str', type=str, default=None)
    args = parser.parse_args()
    inp = args.input
    o = args.output
    k = args.kmer
    q = args.quality
    b = args.baseline



    print(finalize(inp, k, q, b))

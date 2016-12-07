import os
import csv
from collections import defaultdict
import argparse

def get_id_seq_tuple(f):
    while True:
        id = f.readline()
        seq = f.readline()
        if not id: 
            break
        yield (id, seq)

def get_sequence_frequency_for_fasta(path, out_path):
    s_counts = []
    
    with open(path) as f:
        seqs = get_id_seq_tuple(f)
        p_id, p_seq = seqs.next()
        for r in p_seq.strip():
            s_counts.append({
                'primary': r
            })
        for id, seq in seqs:
            for i, r in enumerate(seq.strip()):
                if r != '-':
                    if s_counts[i].get('count') is None:
                        s_counts[i]['count'] = 0
                    s_counts[i]['count'] += 1
                    if r != s_counts[i]['primary']:
                        if s_counts[i].get('variant_count') is None:
                            s_counts[i]['variant_count'] = 0
                        s_counts[i]['variant_count'] += 1
    with open(out_path, 'wb') as out:
        writer = csv.writer(out)
        writer.writerow([p_id])
        writer.writerow(['position', 'residue', 'freq'])
        for i, e in enumerate(s_counts):
            r = e['primary']
            count = e['count']
            variant_count = e.get('variant_count', 0)
            freq = 1.0 - variant_count / float(count)
            writer.writerow([i+1, r, freq])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gathers frequency data from a multifasta file.')
    parser.add_argument('fastapath', metavar='N', nargs=1,
                        help='multifasta imput file')
    parser.add_argument('destpath', nargs=1,
                        help='output path')
    args = parser.parse_args()
    get_sequence_frequency_for_fasta(args.fastapath[0], args.destpath[0])



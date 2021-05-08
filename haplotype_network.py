#!/usr/bin/env python
# coding=utf-8
# __author__ = 'Yunchao Ling'


# import pandas as pd
# import numpy as np
import os
import re
import click
import multiprocessing
from Bio import SeqIO
from tqdm import tqdm, trange

import build.hamming as hamming_cpp

valid_nuc = set(["A", "C", "G", "T"])
pos_freq = []
pos_ref = []


def hamming_py(seq1: str, seq2: str, poss_freq: dict):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    distance = 0
    max_maf = 0
    diff = []
    for i in range(len(seq1)):
        if seq1[i] in valid_nuc and seq2[i] in valid_nuc:
            if seq1[i] != seq2[i]:
                distance += 1
                diff.append(str(i))
                if poss_freq[i] > max_maf:
                    max_maf = poss_freq[i]
        else:
            if seq1[i] != seq2[i]:
                diff.append(str(i))
                # if poss_freq[i] > max_maf:
                #     max_maf = poss_freq[i]
        # if seq1[i] != seq2[i]:
        #     distance += 1
        #     diff.append(str(i))
        #     if poss_freq[i] > max_maf:
        #         max_maf = poss_freq[i]
    if distance == 0:
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                # diff.append(str(i))
                if poss_freq[i] > max_maf:
                    max_maf = poss_freq[i]
    # diff.sort(key=lambda x: x[0])
    # diff = ["%d:%s:%s" % (item[0], item[1], item[2]) for item in diff]
    # diff = "|".join(diff)
    return distance, max_maf, diff


def init_freq(in_file: str):
    global pos_freq
    global pos_ref

    count = 0
    with open(in_file, "r") as infile:
        infile.readline()
        for line in infile:
            splitline = line.rstrip().split("\t")
            pos_ref.append(int(splitline[0]))
            pos_freq.append(float(splitline[2]))
            count += 1


def find_clade(clades: dict, node: int):
    for key in clades:
        if node in clades[key]:
            return key


def seq2geno(seq: str, seqindex: list):
    geno_list = []
    for item in seqindex:
        geno_list.append("%d:%s" % (item[1], seq[item[0]]))
    genos = "|".join(geno_list)
    return genos


def exec_queue_py(iter: int, seqss: list, poss_freq: dict, out_file: str):
    outfile = open(out_file, "a+")
    # for i in tqdm(range(len(seqss)), desc=str(iter)):
    for i in trange(len(seqss), leave=False, desc="H" + str(iter),
                    position=int(multiprocessing.current_process().name.split("-")[1])):
        if iter < i:
            distance, max_maf, diff = hamming_py(seqss[iter], seqss[i], poss_freq)
            outfile.write("%d\t%d\t%d\t%f\t%s\n" % (iter, i, distance, max_maf, ",".join(diff)))
            outfile.flush()
    outfile.close()


@click.command()
@click.argument("in_dir", type=click.Path(exists=True))
def batch_haplotype_network(in_dir: str):
    '''
        批量生成单倍型，并且按照all生成单倍型网络
        输入目录包括pi_pos FASTA文件，freq文件
        输出在同一目录，包括node文件，和net文件
        对所有频率阈值生成node文件，仅对无过滤的all生成net文件
    '''
    for infile in os.listdir(in_dir):
        if infile.startswith("pi_pos"):
            print("正在计算文件%s" % infile)
            pi_pos_file = os.path.join(in_dir, infile)
            freq_file = re.sub(r'pi_pos_(.*?).fasta', r'freq_\1.txt', pi_pos_file)
            if infile.find("all") != -1:
                haplotype_network(pi_pos_file, freq_file, True)
            else:
                haplotype_network(pi_pos_file, freq_file, False)


def haplotype_network(pi_pos_file: str, freq_file: str, draw_net: bool):
    print("生成pattern列表")
    seqs = {}
    seq_count = 0
    seqrecords = SeqIO.parse(pi_pos_file, "fasta")
    for seqrecord in seqrecords:
        seq_count += 1
        # seqs.add(str(seqrecord.seq).upper())
        seq = str(seqrecord.seq).upper()
        if seq in seqs:
            seqs[seq].append(str(seq_count))
        else:
            seqs[seq] = [str(seq_count)]

    init_freq(freq_file)

    print("生成nodes文件")
    seqss = list(seqs)
    seqindex = []
    for i in range(len(seqss[0])):
        seqindex.append([i, pos_ref[i], pos_freq[i]])
    seqindex.sort(key=lambda x: (-x[2], x[1]))
    nodes_file = open(re.sub(r'pi_pos_(.*?).fasta', r'nodes_\1.txt', pi_pos_file), "w")
    genos = {}
    for seq in seqss:
        genos[seq] = seq2geno(seq, seqindex)
        nodes_file.write(genos[seq] + "\t" + "\t".join(seqs[seq]) + "\n")
        nodes_file.flush()
    nodes_file.close()

    if draw_net:
        print("计算海明距离矩阵")
        # df_distance = pd.DataFrame(np.zeros([len(seqss), len(seqss)], dtype=int), index=seqss, columns=seqss)
        # df_max_maf = pd.DataFrame(np.zeros([len(seqss), len(seqss)], dtype=float), index=seqss, columns=seqss)
        # df_diff = pd.DataFrame(np.empty([len(seqss), len(seqss)], dtype=str), index=seqss, columns=seqss)
        # for i in tqdm(range(len(seqss)), desc="line"):
        #     for j in range(len(seqss)):
        #         if i > j:
        #             df_distance.iloc[i, j], df_max_maf.iloc[i, j], df_diff.iloc[i, j] = hamming(seqss[i], seqss[j])
        tempfile = os.path.join(os.path.dirname(pi_pos_file), "candidate_links.txt")
        if os.path.exists(tempfile):
            os.remove(tempfile)
        # p = multiprocessing.Pool(multiprocessing.cpu_count())
        # p = multiprocessing.Pool(multiprocessing.cpu_count(), initializer=tqdm.set_lock,
        #                          initargs=(multiprocessing.RLock(),))
        # for i in range(len(seqss)):
        #     p.apply_async(exec_queue_py, args=(i, seqss, pos_freq, tempfile))
        # p.close()
        # p.join()
        hamming_cpp.exec_queue(seqss, pos_freq, tempfile)

        print("生成候选link列表")
        node_list = []
        # for i in tqdm(range(len(df_distance.index)), desc="line"):
        #     for j in range(len(df_distance.columns)):
        #         if i > j:
        #             # if df_distance.iloc[i, j] > 0:
        #             node_list.append([i + 1, j + 1, df_distance.iloc[i, j], df_max_maf.iloc[i, j], df_diff.iloc[i, j]])
        with open(tempfile, "r") as candi_file:
            for line in tqdm(candi_file, desc="link"):
                line = line.rstrip()
                splitline = line.split("\t")
                node_list.append([int(splitline[0]) + 1, int(splitline[1]) + 1, int(splitline[2]), float(splitline[3]),
                                  splitline[4]])

        print("候选link排序")
        node_list.sort(key=lambda x: (x[2], x[3]))

        print("生成网络")
        net_list = []
        max_length = len(seqss)
        added_nodes = set()
        clades = {}
        current_clade = 0

        for node in tqdm(node_list, desc="link"):
            if len(added_nodes) != max_length:
                if node[0] not in added_nodes:
                    if node[1] not in added_nodes:
                        clades[current_clade] = set([node[0], node[1]])
                        current_clade += 1
                    else:
                        clades[find_clade(clades, node[1])].add(node[0])
                    net_list.append(node)
                else:
                    if node[1] not in added_nodes:
                        clades[find_clade(clades, node[0])].add(node[1])
                        net_list.append(node)
                    else:
                        clade0 = find_clade(clades, node[0])
                        clade1 = find_clade(clades, node[1])
                        if clade0 != clade1:
                            clades[clade0] = clades[clade0].union(clades[clade1])
                            clades.pop(clade1)
                            net_list.append(node)
                added_nodes.add(node[0])
                added_nodes.add(node[1])
            else:
                clade0 = find_clade(clades, node[0])
                clade1 = find_clade(clades, node[1])
                if clade0 != clade1:
                    clades[clade0] = clades[clade0].union(clades[clade1])
                    clades.pop(clade1)
                    net_list.append(node)
                if len(clades) < 2:
                    break

        net_file = open(re.sub(r'pi_pos_(.*?).fasta', r'net_\1.txt', pi_pos_file), "w")
        for link in net_list:
            # net_file.write("%d\t%d\t%d\t%s\n" % (link[0], link[1], link[2], link[4]))
            diff = []
            for i in link[4].split(","):
                i = int(i)
                diff.append([pos_ref[i], seqss[link[0] - 1][i], seqss[link[1] - 1][i]])
            diff.sort(key=lambda x: x[0])
            diff = ["%d:%s:%s" % (item[0], item[1], item[2]) for item in diff]
            diff = "|".join(diff)
            net_file.write("%d\t%d\t%d\t%s\n" % (link[0], link[1], link[2], diff))
            net_file.flush()
        net_file.close()


if __name__ == '__main__':
    batch_haplotype_network()

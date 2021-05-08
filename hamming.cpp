#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <atomic>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;

struct HammingRet
{
	int distance;
	double max_maf;
	vector<string> diff;
};

HammingRet hamming(const string& seq1, const string& seq2, const vector<double>& poss_freq)
{
	//seq1 = seq1.upper()
	//seq2 = seq2.upper()
	int distance = 0;
	double max_maf = 0;
	vector<string> diff;
	for (int i = 0; i < seq1.size(); i++)
	{
		char buf[16];
		if (seq1[i] == 'A' || seq1[i] == 'C' || seq1[i] == 'G' || seq1[i] == 'T')
		{
			if (seq1[i] != seq2[i])
			{
				distance += 1;
				sprintf(buf, "%d", i);
				diff.emplace_back(buf);
				if (poss_freq[i] > max_maf)
					max_maf = poss_freq[i];
			}
		}
		else
		{
			if (seq1[i] != seq2[i])
			{
				sprintf(buf, "%d", i);
				diff.emplace_back(buf);
				if (distance == 0)
				{
					for (int i = 0; i < seq1.size(); i++)
					{
						if (seq1[i] != seq2[i])
							if (poss_freq[i] > max_maf)
								max_maf = poss_freq[i];
					}
				}
			}
		}
	}
	return { distance, max_maf, diff };
}

void exec_queue(const vector<string>& seqss, const vector<double>& poss_freq, const string& out_file)
{
	FILE* out = fopen(out_file.c_str(), "w");
	for (int i = 0; i < seqss.size(); i++)
	{
		for (int j = 0; j < i; j++)
		{
			auto hamming_ret = hamming(seqss[j], seqss[i], poss_freq);
			fprintf(out, "%d\t%d\t%d\t%f\t", j, i, hamming_ret.distance, hamming_ret.max_maf);
			for (int k = 0; k < hamming_ret.diff.size(); k++)
			{
				fprintf(out, "%s%s", hamming_ret.diff[k].c_str(), (k == hamming_ret.diff.size() - 1) ? "\n" : ", ");
			}
		}
	}
}

// int main()
// {
// 	vector<string> seqss(3);
// 	const vector<double> poss_freq(289, 0.5);
// 	seqss[0] = "GGACGCCCGCGTTGCCCATCCCCCCGAGGATCGTCCCGGCGACCGCACCGTCCCCCTACCCGGTGACTCATCCCCCCCGTCCCCCTACCGCGCACCCCGCGACCCCTCCTTTGCCCATAGCCCCCCGTGGTCAAGCCGCGTCCCCCCCATCCTCTCCCCCTGGTCCCACTCGCAGATCTCCGGTTCCGGGTACATCGTCCCCCCCGCACAAGGAACCTGCGCGCGCGCGCAGTTCCCCGGTACGCGCGGCCGGTCTAGCGGCGCGGGGTCTGCGGGGGCACAACGGGC";
// 	seqss[1] = "GGACGCCCGCGTTGCCCATCCCCCCAAGGATCGTCCCGGCGACCGCACCGTCCCCCTACCCGATGACTCATCCCCCCCGCCCCCCTACCGCGCACCCCTCGACCCCTCCTTTGCCCATAGCCCCCCGTGGTCAAGCCGCGTCCCCCCCATCCTCTCCCCCTGGTCCCACTCGCAGATCTCCGGTTCCGGGTACATCGTCCCCCCCGCACAAGGAACCTGCGCGCGCGCGCAGTTCCCCGGTACGCGTGGCCGGTCCAGCGGCGCGGGGTCTGCGGGGGCACAACGTGC";
// 	seqss[2] = "GGACGCCCGCGTTGCCTATCCCCCCAAGGATCGTCCCGGCGACCGCACCGTCCCCCTACCCGGTGACTCATCCCCCCCTCCCCCCTACCGCGCACCCCTCGACCCCTCCTTTGCCCATAGCCCCCCGTGGTCAAGCCGCGTCCCCCCCATCCTCTCCCCCTGGTCCCACTCGCAGATCTCCGGTTCCGGGTACATCGTCCCCCCCGCACAAGGAATCTGCGCGCGCGCGCAGTTCCCCGGTACGCGTGGCCGGTCCAGCGGCGCGGGGTCTGCGGGGGCACAACGTGC";
// 	exec_queue(seqss, poss_freq, "test.txt");
// 	return 0;
// }


PYBIND11_MODULE(hamming, m) {
  m.def("exec_queue", &exec_queue, "Exec queue multi-thread implementation");
}


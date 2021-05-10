#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <atomic>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>
#include <mutex>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdlib.h>
#include <cstdio>

using namespace std;

struct HammingRet
{
	int distance;
	double max_maf;
	vector<int> diff;
};

HammingRet hamming(const string& seq1, const string& seq2, const vector<double>& poss_freq)
{
	int distance = 0;
	double max_maf = 0;
	vector<int> diff;
	diff.reserve(seq1.size());
	for (int i = 0; i < seq1.size(); i++)
	{
		if (seq1[i] == 'A' || seq1[i] == 'C' || seq1[i] == 'G' || seq1[i] == 'T')
		{
			if (seq1[i] != seq2[i])
			{
				distance += 1;
				diff.push_back(i);
				if (poss_freq[i] > max_maf)
					max_maf = poss_freq[i];
			}
		}
		else
		{
			if (seq1[i] != seq2[i])
			{
				diff.push_back(i);
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
	vector<string> seqss_upper(seqss.size());
	for (auto i = 0; i < seqss.size(); i++)
	{
		seqss_upper[i] = seqss[i];
		for (auto j = 0; j < seqss_upper[i].size(); j++)
		{
			char c = seqss_upper[i][j];
			if (c >= 'a' && c <= 'z')
			{
				seqss_upper[i][j] = c - (char)32;
			}
		}
	}
	ofstream out(out_file);
	atomic<size_t> I = 0;
	mutex mu;
	auto fn = [&]()
	{
		size_t n = 6000 / seqss_upper[0].size();
		char* buf = new char[0x200000];
		char* p = buf;
		while (true)
		{
			size_t i_begin = 0;
			do {
				i_begin = I;
			} while (!I.compare_exchange_strong(i_begin, i_begin + n));
			
			if (i_begin >= seqss.size())
				break;
      printf("i_begin %d\n", i_begin);
			auto i_end = i_begin + n;
			if (i_end > seqss_upper.size())
				i_end = seqss_upper.size();
			for (auto j = i_begin + 1; j < seqss_upper.size(); j++)
			{
				for (auto i = i_begin; i < i_end; i++)
				{
					if (i < j)
					{
						auto hamming_ret = hamming(seqss_upper[i], seqss_upper[j], poss_freq);
						p += sprintf(p, "%zu\t%zu\t%d\t%f\t", i, j, hamming_ret.distance, hamming_ret.max_maf);
						if (!hamming_ret.diff.empty())
						{
							for (int k = 0; k < hamming_ret.diff.size(); k++)
							{
                p = _itoa(hamming_ret.diff[k], p, 10);
                while (*p != '\0')
                  p++;
								// p += sprintf(p, "%d", hamming_ret.diff[k]);
								if (k != hamming_ret.diff.size() - 1)
									*p++ = ',';
							}
						}
						*p++ = '\n';
						if (p - buf > 0x1f0000)
						{
							*p = '\0';
							mu.lock();
							out << buf;
							mu.unlock();
							p = buf;
						}
					}
				}
			}
		}
		if (p != buf)
		{
			*p = '\0';
			mu.lock();
			out << buf;
			mu.unlock();
			p = buf;
		}
		delete[] buf;
	};
	vector<thread> threads(8); // 线程数这里改
	for (int i = 0; i < threads.size(); i++)
		threads[i] = thread(fn);
	for (int i = 0; i < threads.size(); i++)
		threads[i].join();
	out.close();
}

// int main()
// {
// 	int n = 2000;
// 	double rate = 0.8;
// 	uniform_int_distribution<int> u1(160000, 999999);
// 	uniform_real_distribution<double> u2(0.0, 1.0);
// 	uniform_int_distribution<int> u3(0, 288);
// 	uniform_int_distribution<int> u4(0, 3);
// 	default_random_engine e1(time(0) + clock());
// 	default_random_engine e2(time(0) + clock());
// 	default_random_engine e3(time(0) + clock());
// 	default_random_engine e4(time(0) + clock());
// 	map<string, string> M;

// 	char key[64];
// 	sprintf(key, ">OEAV%06d", 2);
// 	string base_val = "aGACGCCCGCGTTGCCCATCCCCCCGAGGATCGTCCCGGCGACCGCACCGTCCCCCTACCCGGTGACTCATCCCCCCCGTCCCCCTACCGCGCACCCCGCGACCCCTCCTTTGCCCATAGCCCCCCGTGGTCAAGCCGCGTCCCCCCCATCCTCTCCCCCTGGTCCCACTCGCAGATCTCCGGTTCCGGGTACATCGTCCCCCCCGCACAAGGAACCTGCGCGCGCGCGCAGTTCCCCGGTACGCGCGGCCGGTCTAGCGGCGCGGGGTCTGCGGGGGCACAACGGGC";
// 	M.insert(make_pair(key, base_val));
// 	while (M.size() < n)
// 	{
// 		sprintf(key, ">OEAV%06d", u1(e1));
// 		string val = base_val;
// 		while (u2(e2) < rate)
// 		{
// 			static const char AGCT[4] = { 'A', 'G', 'C', 'T' };
// 			val[u3(e3)] = AGCT[u4(e4)];
// 		}
// 		M.insert(make_pair(key, val));
// 	}
// 	vector<string> seqss;
// 	for (auto it = M.begin(); it != M.end(); ++it)
// 	{
// 		seqss.push_back(it->second);
// 	}
// 	const vector<double> poss_freq(289, 0.5);
// 	cout << "start" << endl;
// 	clock_t t0 = clock();
// 	exec_queue(seqss, poss_freq, "test.txt");
// 	cout << clock() - t0 << endl;
// 	return 0;
// }

PYBIND11_MODULE(hamming, m) {
  m.def("exec_queue", &exec_queue, "Exec queue multi-thread implementation");
}

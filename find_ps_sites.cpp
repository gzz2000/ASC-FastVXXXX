#include <vector>
#include <utility>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include <pybind11/stl_bind.h>

// PYBIND11_MAKE_OPAQUE(std::tuple<std::vector<int>, std::vector<int>>);
// PYBIND11_MAKE_OPAQUE(std::vector<int>);
// PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
// PYBIND11_MAKE_OPAQUE(std::string);

const int num_nc = 8;

inline int idx_nc(int c) {
  switch(c) {
  case 'a': return 0;
  case 't': return 1;
  case 'g': return 2;
  case 'c': return 3;
  case 'A': return 4;
  case 'T': return 5;
  case 'C': return 6;
  case 'G': return 7;
  default: return -1;
  }
}

std::tuple<std::vector<int>, std::vector<int>> find_ps_sites(const std::vector<std::string> &alignment, const std::string &seq, double remove_fq, int epis_base)
{
  // __builtin_debugtrap();
  // epis_base += 1;
  std::vector<int> uniq_nc(alignment[0].length() * num_nc, 0);
  for(int i = 0; i < alignment.size(); ++i) {
    for(int j = 0; j < alignment[0].length(); ++j) {
      int t = idx_nc(alignment[i][j]);
      if(t < 0) continue;
      ++uniq_nc[j * num_nc + t];
    }
  }
  // epis_base -= 1;
  std::vector<int> all_pis, rm_pis;
  for(int j = 0; j < alignment[0].length(); ++j) {
    int effective = 0, cur = 0, mx = -1, mxk = -1, sum = 0;
    for(int k = 0; k < num_nc; ++k) {
      int v = uniq_nc[j * num_nc + k];
      if(v) ++effective;
      if(v > 1) ++cur;
      if(v > mx) { mx = v; mxk = k; }
      sum += v;
    }
    if(effective <= 1 || sum <= epis_base) continue;
    if(cur > 1 && seq[j] != '-') {
      all_pis.push_back(j);
      double pos_freq = 1. - (double)mx / sum;
      if(pos_freq > remove_fq) rm_pis.push_back(j);
    }
  }
  return std::make_tuple(std::move(all_pis), std::move(rm_pis));
}

PYBIND11_MODULE(find_ps_sites, m) {
  // py::bind_vector<std::vector<std::string>>(m, "VectorString");
  m.def("find_ps_sites", &find_ps_sites, "Step 1 find_ps_sites CPP implementation");
}

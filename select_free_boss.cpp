#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <algorithm>
#include "stdlib_printing.hh"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::cout;
using std::endl;

class SelectFreeBOSS{
public:
  SelectFreeBOSS(const vector<string>& input, int k);
  // One bit vector for each character. Can be empty bit vector if the character does not occur
  string SBWT[256];
  // SBWT['A'], SBWT['C'], SBWT['G'], SBWT['T'] are all bit vectors
  vector<int> C; // C-array (cumulative character counts)
  int n_nodes;
};

// true if S is colexicographically-smaller than T
bool colex_compare(const string& S, const string& T) {
  for (int i_s = S.size() - 1, i_t = T.size() - 1;; --i_s, --i_t){
    // One of the strings is a suffix of the other. Return the shorter.
    if(i_s < 0 || i_t < 0) return S.size() < T.size();
    if(S[i_s] != T[i_t]) return S[i_s] < T[i_t];
  }
}

// Counts the number of occurrence of symbol in array[0..position)
int64_t Rank(const string& S, char symbol, int64_t position){
  int64_t ans = 0;
  for(int64_t i = 0; i < position; ++i)
    if(S[i] == symbol) ++ans;
  return ans;
}

// Returns position i such that array[i] == symbol and
// symbol occurs `count` times in array[0..i]
// Using capital S in the name because select conflicts with the standard library
int64_t Select(const string& S, char symbol, int64_t count){
  int64_t occurrences_seen = 0;
  for(int64_t i = 0; i < S.size(); ++i){
    if(S[i] == symbol) occurrences_seen++;
    if(occurrences_seen == count) return i;
  }
  throw std::range_error("Select reached end of string before finding the requested number of symbols");
}

vector<int> cumulative_sum_of_counts(const vector<int>& counts) {
  vector<int> C(counts);
  for (int i = 1; i < C.size(); ++i)
    C[i] += C[i-1];
  return C;
}

template <typename T>
inline void shift_vector_to_the_right_by_1(vector<T>& v){
  for(int i = v.size() - 1; i > 0; --i)
    v[i] = v[i-1];
  v[0] = 0;
}

vector<int> construct_C(const vector<int>& counts){
  vector<int> C = cumulative_sum_of_counts(counts);
  shift_vector_to_the_right_by_1(C); // we shift to follow the definition
  return C;
}

// Edge-centric definition.
// k is the length of node labels.
SelectFreeBOSS::SelectFreeBOSS(const vector<string>& input, int k){
  map<string, std::pair<set<char>, set<char>>, decltype(colex_compare)*> kmers(colex_compare); // k-mer -> (incoming labels, outgoing labels)

  // TODO: ensure that root node exists
  // TODO: avoid adding redundant dummies

  for(const string& S : input){
    assert(S.size() >= k);

    // Dummies
    for(int i = 0; i <= k; i++){
      string prefix = string(k-i, '$') + S.substr(0,i);
      if(i != 0) kmers[prefix].first.insert('$'); // In
      if(i < k) kmers[prefix].second.insert(S[i]); // Out
    }

    // Non-dummies
    for(int i = 0; i < S.size()-k+1; i++){
      string kmer = S.substr(i,k);
      if(i > 0) kmers[kmer].first.insert(S[i-1]); // In
      if(i + k < S.size()) kmers[kmer].second.insert(S[i+k]); // Out
    }
  }

  this->n_nodes = kmers.size();

  // Add dollars
  for(auto& keyval : kmers) {
    if(keyval.second.second.size() == 0){
      keyval.second.second.insert('$'); // Outgoing dollar
      kmers.begin()->second.first.insert('$'); // Incoming dollar to root
    }
  }

  string SBWT[256];
  SBWT['A'].resize(this->n_nodes, '0');
  SBWT['C'].resize(this->n_nodes, '0');
  SBWT['G'].resize(this->n_nodes, '0');
  SBWT['T'].resize(this->n_nodes, '0');
  string F_column;
  int kmer_index = 0;
  for(auto& keyval : kmers) {
    cout << keyval.first << " " << keyval.second << endl;
    for(char c : keyval.second.second){
      SBWT[c][kmer_index] = '1';
      F_column += c;
    }
    kmer_index++;
  }

  this->SBWT['A'] = SBWT['A'];
  this->SBWT['C'] = SBWT['C'];
  this->SBWT['G'] = SBWT['G'];
  this->SBWT['T'] = SBWT['T'];

  // Add minus marks to SBWT
  std::sort(F_column.begin(), F_column.end());
  map<char,int> labels_seen;
  int F_index = 0;
  vector<int> counts(256);
  for(char c : F_column) counts[c]++; // Will subtract minus-characters below
  for(auto& keyval : kmers) {
    int indegree = keyval.second.first.size();
    char c = F_column[F_index];
    if(c != '$'){ // Add minuses, but not for dollars
      for(int i = 1; i < indegree; i++){ // All but the first in-edge
        this->SBWT[c][Select(SBWT[c], '1', labels_seen[c] + i + 1)] = '0'; // Turn off the bit
        counts[c]--;
      }
    }
    labels_seen[c] += indegree;
    F_index += indegree;
  }

  this->C = construct_C(counts);

  cout << "SBWT[\'A\'] = " << this->SBWT['A'] << endl;
  cout << "SBWT[\'C\'] = " << this->SBWT['C'] << endl;
  cout << "SBWT[\'G\'] = " << this->SBWT['G'] << endl;
  cout << "SBWT[\'T\'] = " << this->SBWT['T'] << endl;
  cout << this->C << endl;
}


int search(SelectFreeBOSS& boss, const string& kmer){
  int node_left = 0;
  int node_right = boss.n_nodes-1;
  for(int i = 0; i < kmer.size(); i++){
    char c = kmer[i];
    node_left = boss.C[c] + Rank(boss.SBWT[c], '1', node_left);
    node_right = boss.C[c] + Rank(boss.SBWT[c], '1', node_right+1) - 1;
    if(node_left > node_right) return -1; // Not found
  }
  assert(node_left == node_right);
  return node_left;
}


int main(){
  vector<string> input = {"GAAGCCGCCATTCCATAGTGAGTCCTTCGTCTGTGACTATCTGTGCCAGATCGTCTAGCAAACTGCTGATCCAGTTTATCTCACCAAATTATAGCCGTACAGACCGAAATCTTAAGTCATATCACGCGACTAGGCTCAGCTTTATTTTTGTGGTCATGGGTTTTGGTCCGCCCGAGCGGTGCAGCCGATTAGGACCATGT"};
  int k = 4;
  SelectFreeBOSS boss(input, k);
  set<string, decltype(colex_compare)*> kmers(colex_compare);
  for(string& S : input)
    for(int i = 0; i < S.size()-k+1; i++)
      kmers.insert(S.substr(i,k));


  for(string kmer : kmers){
    int result = search(boss, kmer);
    cout << result << endl;
  }
}

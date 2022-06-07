#include <iostream>
#include <map>
#include <string>
#include <cassert>
#include <algorithm>
#include "stdlib_printing.hh"

using namespace std;

struct SelectFreeBOSS{
  string GBWT[256]; // One bit vector for each character. Can be empty bit vector if the character does not occur
  //string GBWT_A; // Bit vector
  //string GBWT_C; // Bit vector
  //string GBWT_G; // Bit vector
  //string GBWT_T; // Bit vector
  vector<int> C; // C-array (cumulative character counts)
  int n_nodes;
};

struct colex_compare {
  // true if S is colexicographically-smaller than T
  bool operator()(const std::string& S, const std::string& T) const {
    int i = 0;
    while(true){
      if(i == S.size() || i == T.size()){
        // One of the strings is a suffix of the other. Return the shorter.
        if(S.size() < T.size()) return true;
        else return false;
      }
      if(S[S.size()-1-i] < T[T.size()-1-i]) return true;
      if(S[S.size()-1-i] > T[T.size()-1-i]) return false;
      i++;
    }
  }
};

// Counts the number of occurrence of symbol in array[0..position)
int64_t Rank(const string& S, char symbol, int64_t position){
  int64_t ans = 0;
  for(int64_t i = 0; i < position; i++)
    if(S[i] == symbol) ans++;
  return ans;
}

// Returns position i such that array[i] == symbol and
// symbol occurs `count` times in array[0..i]
// Using capital S in the name because select conflicts with the standard library
int64_t Select(const string& S, char symbol, int64_t count){
  int64_t occurrences_seen = 0;
  for(int64_t i = 0; ; i++){
    if(S[i] == symbol) occurrences_seen++;
    if(occurrences_seen == count) return i;
  }
}

vector<int> char_counts_to_C_array(const vector<int>& counts){
  vector<int> C(256); // Cumulative sum of counts

  // Compute cumulative sum of counts
  for(int i = 0; i < (int)C.size(); i++){
    C[i] = counts[i];
    if(i > 0) C[i] += C[i-1];
  }

  // Shift C to the right by one because that's how it's defined
  for(int i = 256-1; i >= 0; i--){
    if(i == 0) C[i] = 0;
    else C[i] = C[i-1];
  }

  return C;
}


// Edge-centric definition.
// k is the length of node labels.
SelectFreeBOSS construct(const vector<string>& input, int k){
  map<string, pair<set<char>, set<char>>, colex_compare> kmers; // k-mer -> (incoming labels, outgoing labels)

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

  SelectFreeBOSS boss;
  boss.n_nodes = kmers.size();

  // Add dollars
  for(auto& keyval : kmers) {
    if(keyval.second.second.size() == 0){
      keyval.second.second.insert('$'); // Outgoing dollar
      kmers.begin()->second.first.insert('$'); // Incoming dollar to root
    }
  }

  string GBWT[256];
  GBWT['A'].resize(boss.n_nodes, '0');
  GBWT['C'].resize(boss.n_nodes, '0');
  GBWT['G'].resize(boss.n_nodes, '0');
  GBWT['T'].resize(boss.n_nodes, '0');
  string F_column;
  int kmer_index = 0;
  for(auto& keyval : kmers) {
    cout << keyval.first << " " << keyval.second << endl;
    for(char c : keyval.second.second){
      GBWT[c][kmer_index] = '1';
      F_column += c;
    }
    kmer_index++;
  }

  boss.GBWT['A'] = GBWT['A'];
  boss.GBWT['C'] = GBWT['C'];
  boss.GBWT['G'] = GBWT['G'];
  boss.GBWT['T'] = GBWT['T'];

  // Add minus marks to GBWT
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
        boss.GBWT[c][Select(GBWT[c], '1', labels_seen[c] + i + 1)] = '0'; // Turn off the bit
        counts[c]--;
      }
    }
    labels_seen[c] += indegree;
    F_index += indegree;
  }

  boss.C = char_counts_to_C_array(counts);

  cout << "GBWT[\'A\'] = " << boss.GBWT['A'] << endl;
  cout << "GBWT[\'C\'] = " << boss.GBWT['C'] << endl;
  cout << "GBWT[\'G\'] = " << boss.GBWT['G'] << endl;
  cout << "GBWT[\'T\'] = " << boss.GBWT['T'] << endl;
  cout << boss.C << endl;
  return boss;
}


int search(SelectFreeBOSS& boss, const string& kmer){
  int node_left = 0;
  int node_right = boss.n_nodes-1;
  for(int i = 0; i < kmer.size(); i++){
    char c = kmer[i];
    node_left = boss.C[c] + Rank(boss.GBWT[c], '1', node_left);
    node_right = boss.C[c] + Rank(boss.GBWT[c], '1', node_right+1) - 1;
    if(node_left > node_right) return -1; // Not found
  }
  assert(node_left == node_right);
  return node_left;
}


int main(){
  vector<string> input = {"GAAGCCGCCATTCCATAGTGAGTCCTTCGTCTGTGACTATCTGTGCCAGATCGTCTAGCAAACTGCTGATCCAGTTTATCTCACCAAATTATAGCCGTACAGACCGAAATCTTAAGTCATATCACGCGACTAGGCTCAGCTTTATTTTTGTGGTCATGGGTTTTGGTCCGCCCGAGCGGTGCAGCCGATTAGGACCATGT"};
  int k = 4;
  SelectFreeBOSS boss = construct(input, k);
  set<string, colex_compare> kmers;
  for(string& S : input)
    for(int i = 0; i < S.size()-k+1; i++)
      kmers.insert(S.substr(i,k));


  for(string kmer : kmers){
    int result = search(boss, kmer);
    cout << result << endl;
  }
}

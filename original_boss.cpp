#include <iostream>
#include <map>
#include <string>
#include <cassert>
#include <algorithm>
#include "stdlib_printing.hh"

using namespace std;

struct BOSS{

    string GBWT; // Generalized BWT
    string LAST; // Bit vector 'last'
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
BOSS construct(const vector<string>& input, int k){
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

    BOSS boss;
    boss.n_nodes = kmers.size();

    // Add dollars
    for(auto& keyval : kmers) {
        if(keyval.second.second.size() == 0){
            keyval.second.second.insert('$'); // Outgoing dollar
            kmers.begin()->second.first.insert('$'); // Incoming dollar to root
        }
    }

    // Construct GBWT with only capital letters, and LAST
    for(auto& keyval : kmers) {
        cout << keyval.first << " " << keyval.second << endl;
        for(char c : keyval.second.second) boss.GBWT += c;
        while(boss.LAST.size()+1 < boss.GBWT.size()) boss.LAST += '0';
        boss.LAST += '1';
    }

    // Add minus marks to GBWT
    string F_column = boss.GBWT;
    std::sort(F_column.begin(), F_column.end());
    vector<bool> minus_marks(boss.GBWT.size());
    map<char,int> labels_seen;
    int F_index = 0;
    for(auto& keyval : kmers) {
        int indegree = keyval.second.first.size();        
        char c = F_column[F_index];
        if(c != '$'){ // Add minuses, but not for dollars
            for(int i = 1; i < indegree; i++){ // All but the first in-edge
                minus_marks[Select(boss.GBWT, c, labels_seen[c] + i + 1)] = 1; // Mark the corresponding out-edge
            }
        }
        labels_seen[c] += indegree;
        F_index += indegree;
    }

    // Apply minus marks
    for(int i = 0; i < boss.GBWT.size(); i++)
        if(minus_marks[i]) boss.GBWT[i] = tolower(boss.GBWT[i]);

    vector<int> counts(256);
    for(char c : boss.GBWT) counts[c]++;
    boss.C = char_counts_to_C_array(counts);

    cout << boss.GBWT << endl;
    cout << boss.C << endl;    
    cout << boss.LAST << endl;
    return boss;

}

int search(BOSS& boss, const string& kmer){
    int node_left = 0;
    int node_right = boss.n_nodes-1;
    for(int i = 0; i < kmer.size(); i++){
        int GBWT_left = node_left == 0 ? 0 : Select(boss.LAST, '1', node_left) + 1; // End of previous node +1.
        int GBWT_right = Select(boss.LAST, '1', node_right+1);
        char c = kmer[i];
        node_left = boss.C[c] + Rank(boss.GBWT, c, GBWT_left);
        node_right = boss.C[c] + Rank(boss.GBWT, c, GBWT_right+1) - 1;
        if(node_left > node_right) return -1; // Not found
    }
    assert(node_left == node_right);
    return node_left;
}

int main(){
    vector<string> input = {"TACGACGTCGACT"};
    BOSS boss = construct(input, 3);
    vector<string> kmers = {"CGA", "GAC", "TAC", "GTC", "ACG", "TCG", "ACT", "CGT"};
    vector<int> colex_ranks = {1,3,4,5,6,7,9,10};
    for(int i = 0; i < kmers.size(); i++){
        int result = search(boss, kmers[i]);
        if(result != colex_ranks[i]){
            cerr << "Error: k-mer " << kmers[i] << " gave wrong answer (" << result << " vs " << colex_ranks[i] << endl;
        } else{
            cout << colex_ranks[i] << " " << result << endl;
        }
    }
}
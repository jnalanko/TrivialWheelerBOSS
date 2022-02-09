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


// Edge-centric definition.
// k is the length of node labels.
BOSS construct(const vector<string>& input, int k){
    map<string, pair<set<char>, set<char>>, colex_compare> kmers; // k-mer -> (incoming labels, outgoing labels)
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

    cout << boss.GBWT << endl;
    cout << boss.LAST << endl;
    return boss;

}

int main(){
    vector<string> input = {"TACGACGTCGACT"};
    construct(input, 3);
}
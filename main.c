#include "inttypes.h"
#include "assert.h"
#include "stdlib.h"
#include "stdio.h"

typedef struct WheelerBOSS{

    char* I; // String of digits '0' and '1'
    char* O; // String of digits '0' and '1'
    char* GBWT; // Generalized BWT = string characters 'A', 'C', 'G' and 'T'
    int64_t* C; // Has constant length 256
    int64_t n_nodes;
    int64_t n_edges;

} WheelerBOSS;

// Counts the number of occurrence of symbol in array[0..position)
int64_t Rank(char* array, char symbol, int64_t position){
    int64_t ans = 0;
    for(int64_t i = 0; i < position; i++)
        if(array[i] == symbol) ans++;
    return ans;
}

// Returns position i such that array[i] == symbol and
// symbol occurs `count` times in array[0..i]
// Using capital S in the name because select conflicts with the standard library
int64_t Select(char* array, char symbol, int64_t count){
    int64_t occurrences_seen = 0;
    for(int64_t i = 0; ; i++){
        if(array[i] == symbol) occurrences_seen++;
        if(occurrences_seen == count) return i;
    }
}

int64_t search(WheelerBOSS* boss, const char* kmer, int64_t k){
    int64_t left = 0;
    int64_t right = boss->n_nodes-1;
    for(int64_t i = 0; i < k; i++){
        char c = kmer[i];

        int64_t start = Select(boss->O, '1', left+1) - left;

        int64_t end;
        if(right == boss->n_nodes-1) end = boss->n_edges-1; // Last position
        else end = Select(boss->O, '1', right+2) - right - 2;

        if(end < start) return -1; // K-mer not found

        int64_t edge_left = Rank(boss->GBWT, c, start);
        int64_t edge_right = Rank(boss->GBWT, c, end+1);

        if(edge_left == edge_right) return -1; // K-mer not found

        int64_t edge_wheeler_left = boss->C[c] + edge_left;
        int64_t edge_wheeler_right = boss->C[c] + edge_right - 1;

        left = Rank(boss->I, '1', Select(boss->I, '0', edge_wheeler_left+1))-1;
        right = Rank(boss->I, '1', Select(boss->I, '0', edge_wheeler_right+1))-1;

    }
    
    assert(left == right); // If this is wrong the WheelerBOSS is corrupt or the k is wrong.
    return left; // The colexicographic rank of the k-mer
}

int main(int argc, char** argv){
    // Construct the example from 
    // Alanko, J., et al. "Buffering Updates Enables Efficient Dynamic de Bruijn Graphs." (2021).
    char* I = "11010101001010101010101010";
    char* O = "10100101110101010101001010";
    char* GBWT = "ACGCAGGTTACAA";
    int64_t* C = calloc(256, sizeof(int64_t));
    C['A'] = 0;
    C['C'] = 5;
    C['G'] = 8;
    C['T'] = 11;
    int64_t n_nodes = 13; // Includes technical dummy nodes
    int64_t n_edges = 13; // Includes technical dummy edges. Happens to be the same as n_nodes in this example
    WheelerBOSS boss = {I, O, GBWT, C, n_nodes, n_edges};

    // Test that all present k-mers are found

    // Prepare the test
    int64_t test_kmer_count = 9;
    const char *kmers[] = {"ACA", "CGA", "GTA", "CAC", "CGC", "ACG", "GCG", "AGT", "CGT"};
    const int64_t kmer_colex_ranks[] = {2,3,4,6,7,9,10,11,12};

    // Run the test
    for(int64_t i = 0; i < test_kmer_count; i++){
        int64_t colex = search(&boss, kmers[i], 3);
        printf("%s: %" PRId64 "\n", kmers[i], colex);
        if(colex != kmer_colex_ranks[i]){
            printf("ERROR: Query returned wrong answer for kmer %s", kmers[i]);
        }
    }

    // Try to search for a k-mer that is not in the structure
    if(search(&boss, "TGA", 3) != -1){
        printf("ERROR: query found k-mer TGA even though it's not supposed to");
    }
}

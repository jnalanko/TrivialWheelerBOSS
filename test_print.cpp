/*
Tests the functions inside of stdlib_printing.hh
*/


#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <iostream>
#include "stdlib_printing.hh"


using namespace std;

int main() {
  unordered_map<int, string> um = {
    {2, "two"},
    {3, "three"},
    {5, "five"}
  };
  cout << um << endl;

  map<int, string> m = {
    {2, "two"},
    {3, "three"},
    {5, "five"}
  };
  cout << m << endl;

  vector<int> v = {2, 2, 3, 5};
  cout << v << endl;

  set<int> s = {2, 2, 3, 5};
  cout << s << endl;

  multiset<int> ms = {2, 2, 3, 5};
  cout << ms << endl;

  pair<int, string> p = {2, "two"};
  cout << p << endl;
}

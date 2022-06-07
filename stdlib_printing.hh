#pragma once

#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <iostream>

using std::set;
using std::map;
using std::vector;
using std::unordered_map;
using std::ostream;
using std::pair;
using std::multiset;
using std::string;


template <typename Structure, typename Function>
ostream& print_all(
  ostream& os,
  const string prefix,
  const Structure& v,
  const Function print_item,
  const string suffix
) {
  os << prefix;
  for(auto it = v.begin(); it != v.end(); it++) {
    if(it != v.begin()) os << ", ";
    print_item(os, *it);
  }
  os << suffix;
  return os;
}

template <typename T>
void print_individual(ostream& os, const T& item) {
  os << item;
}

template <typename T>
void print_first_second(ostream& os, const T& item) {
  os << item.first << ": " << item.second;
}


template <typename S, typename T>
ostream& operator<<(ostream& os, const unordered_map<S,T>& v){
  return print_all<typeof(v)>(os, "[", v, print_first_second<pair<S, T>>, "]");
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const map<S,T>& v){
  return print_all<typeof(v)>(os, "{", v, print_first_second<pair<S, T>>, "}");
}

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v){
  return print_all<typeof(v)>(os, "[", v, print_individual<T>, "]");
}

template <typename T>
ostream& operator<<(ostream& os, const set<T>& v){
  return print_all<typeof(v)>(os, "[", v, print_individual<T>, "]");
}

template <typename T>
ostream& operator<<(ostream& os, const multiset<T>& v){
  return print_all<typeof(v)>(os, "[", v, print_individual<T>, "]");
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const pair<S,T>& x){
  os << "(" << x.first << ", " << x.second << ")";
  return os;
}

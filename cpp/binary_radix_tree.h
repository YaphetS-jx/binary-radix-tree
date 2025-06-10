#ifndef BINARY_RADIX_TREE_H
#define BINARY_RADIX_TREE_H
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <array>
#include "morton_code.h"

inline int count_leading_zeros(std::uint64_t x, int num_bits){
    if(x==0) return num_bits;
    int bit_len = 64 - __builtin_clzll(x);
    return num_bits - bit_len;
}

inline int delta(const std::vector<std::uint64_t>& morton_codes, int i, int j, int num_bits){
    int n = morton_codes.size();
    if(j<0 || j>=n) return -1;
    return count_leading_zeros(morton_codes[i]^morton_codes[j], num_bits);
}

inline std::pair<int,int> determine_range(const std::vector<std::uint64_t>& morton_codes, int i, int num_bits){
    int n=morton_codes.size();
    int d = (delta(morton_codes,i,i+1,num_bits)-delta(morton_codes,i,i-1,num_bits))>=0?1:-1;
    int delta_min = delta(morton_codes,i,i-d,num_bits);
    int lmax=2;
    while(delta(morton_codes,i,i+lmax*d,num_bits)>delta_min) lmax*=2;
    int l=0; int t=lmax/2;
    while(t>=1){ if(delta(morton_codes,i,i+(l+t)*d,num_bits)>delta_min) l+=t; t/=2; }
    int j=i+l*d; return {std::min(i,j),std::max(i,j)};
}

inline int find_split(const std::vector<std::uint64_t>& morton_codes, int delta_node, int first, int last, int num_bits){
    int split=first; int stride=last-first;
    while(stride>1){
        stride=(stride+1)/2;
        int mid=split+stride;
        if(mid<last && delta(morton_codes,first,mid,num_bits)>delta_node) split=mid;
    }
    return split;
}

inline std::pair<std::uint64_t,std::uint64_t> find_range_3d(std::uint64_t prefix, int num_bits_prefix, int num_bits){
    std::uint64_t mask1 = ((1ULL<<num_bits_prefix)-1)<<(num_bits-num_bits_prefix);
    std::uint64_t mask2 = mask1 ^ ((1ULL<<num_bits)-1);
    std::uint64_t lower_bound = prefix & mask1;
    std::uint64_t upper_bound = prefix | mask2;
    return {lower_bound,upper_bound};
}

struct Node{ int left; int right; std::uint64_t lower_bound; std::uint64_t upper_bound; };

class BinaryRadixTree{
public:
    BinaryRadixTree(const std::vector<std::uint64_t>& morton_codes, int num_bits):morton_codes_(morton_codes),num_bits_(num_bits){ nodes_=build_tree(); }
    template <class Check>
    std::vector<int> traverse_tree(const Check& check) const {
        std::vector<int> results; int n=morton_codes_.size(); std::vector<int> stack; stack.push_back(n); while(!stack.empty()){ int node_idx=stack.back(); stack.pop_back(); if(node_idx<n){ if(check.check_inside(morton_codes_[node_idx])) results.push_back(node_idx); } else { const Node& node=nodes_[node_idx-n]; if(check.check_intersect(node.lower_bound,node.upper_bound)){ stack.push_back(node.left); stack.push_back(node.right);} }} std::sort(results.begin(),results.end()); return results; }
    template <class Check>
    std::unordered_map<int,std::vector<std::array<int,3>>> traverse_tree_distance(const Check& check) const {
        std::unordered_map<int,std::vector<std::array<int,3>>> results; int n=morton_codes_.size(); std::vector<int> stack; stack.push_back(n); while(!stack.empty()){ int node_idx=stack.back(); stack.pop_back(); if(node_idx<n){ auto temp=check.check_inside(morton_codes_[node_idx]); if(!temp.empty()) results[node_idx]=temp; } else { const Node& node=nodes_[node_idx-n]; if(check.check_intersect(node.lower_bound,node.upper_bound)){ stack.push_back(node.left); stack.push_back(node.right); } } } return results; }
private:
    std::vector<Node> build_tree(){ int n=morton_codes_.size(); std::vector<Node> nodes(n-1); for(int i=0;i<n-1;++i){ auto [first,last]=determine_range(morton_codes_,i,num_bits_); int delta_node=delta(morton_codes_,first,last,num_bits_); int split=find_split(morton_codes_,delta_node,first,last,num_bits_); auto range=find_range_3d(morton_codes_[split],delta_node,num_bits_); int left=split==first?split:split+n; int right=split+1==last?split+1:split+1+n; nodes[i]={left,right,range.first,range.second}; } return nodes; }
    std::vector<std::uint64_t> morton_codes_;
    int num_bits_;
    std::vector<Node> nodes_;
};

#endif

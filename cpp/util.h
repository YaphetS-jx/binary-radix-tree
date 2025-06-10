#ifndef UTIL_H
#define UTIL_H
#include <vector>

inline std::vector<int> prefix_sum(const std::vector<int>& arr){
    std::vector<int> pref(arr.size()+1,0);
    for(size_t i=0;i<arr.size();++i) pref[i+1]=pref[i]+arr[i];
    return pref;
}

#endif

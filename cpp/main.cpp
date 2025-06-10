#include "binary_radix_tree.h"
#include "check.h"
#include <iostream>

int main(){
    // Example usage with simple points
    std::vector<std::uint64_t> morton_codes = {0,1,2,3,4,5,6,7};
    BinaryRadixTree tree(morton_codes, 30);
    simple_check check(0,7);
    auto res = tree.traverse_tree(check);
    for(int idx: res) std::cout << idx << " ";
    std::cout << std::endl;
    return 0;
}

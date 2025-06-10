#ifndef MORTON_CODE_H
#define MORTON_CODE_H
#include <cstdint>
#include <tuple>
#include <iostream>

inline void print_bits(std::uint64_t num, int num_bits) {
    std::string bits;
    for(int i=num_bits-1;i>=0;--i) bits += ((num>>i)&1)?'1':'0';
    std::cout << bits << std::endl;
}

inline std::uint64_t interleave_bits(std::uint64_t x, std::uint64_t y, std::uint64_t z, int num_bits_per_dim) {
    std::uint64_t result = 0;
    for(int i=0;i<num_bits_per_dim;++i) {
        result |= ((x & (1ULL<<i)) << (2*i)) |
                  ((y & (1ULL<<i)) << (2*i+1)) |
                  ((z & (1ULL<<i)) << (2*i+2));
    }
    return result;
}

inline std::tuple<std::uint64_t,std::uint64_t,std::uint64_t> deinterleave_bits(std::uint64_t morton_code, int num_bits_per_dim) {
    std::uint64_t x=0,y=0,z=0;
    for(int i=0;i<num_bits_per_dim;++i){
        x |= (morton_code & (1ULL<<(3*i))) >> (2*i);
        y |= (morton_code & (1ULL<<(3*i+1))) >> (2*i+1);
        z |= (morton_code & (1ULL<<(3*i+2))) >> (2*i+2);
    }
    return {x,y,z};
}

inline std::uint64_t fractional_to_binary(double num, int num_bits) {
    std::uint64_t result=0; int n=0; while(num>0 && n<num_bits){ num *=2; if(num>=1.0){ result|=1ULL<<(num_bits-1-n); num -=1.0; } n++; }
    return result;
}

inline std::uint64_t coordinate_to_morton_code(double num, int num_bits){ return fractional_to_binary(num, num_bits-1); }

inline double binary_to_fractional(std::uint64_t binary, int num_bits){ double result=0.0; for(int i=0;i<num_bits;++i){ result += (binary & 1ULL) * 1.0 / (1ULL<<(num_bits-i)); binary >>=1; } return result; }

inline double morton_code_to_coordinate(std::uint64_t morton_code, int num_bits){ return binary_to_fractional(morton_code, num_bits-1); }

inline std::uint64_t coordinate_to_morton_code_3d(double x,double y,double z,int num_bits_per_dim){
    std::uint64_t x_mc = coordinate_to_morton_code(x,num_bits_per_dim);
    std::uint64_t y_mc = coordinate_to_morton_code(y,num_bits_per_dim);
    std::uint64_t z_mc = coordinate_to_morton_code(z,num_bits_per_dim);
    return interleave_bits(x_mc,y_mc,z_mc,num_bits_per_dim);
}

#endif

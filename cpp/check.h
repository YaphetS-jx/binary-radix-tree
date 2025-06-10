#ifndef CHECK_H
#define CHECK_H
#include <cstdint>
#include <vector>
#include <cmath>
#include <array>
#include "morton_code.h"

inline bool is_less_than_or_equal(std::uint64_t morton1, std::uint64_t morton2){
    const std::uint64_t x_mask=0x9249249249249249ULL;
    const std::uint64_t y_mask=0x2492492492492492ULL;
    const std::uint64_t z_mask=0x4924924924924924ULL;
    std::uint64_t x1=morton1 & x_mask; std::uint64_t x2=morton2 & x_mask;
    std::uint64_t y1=morton1 & y_mask; std::uint64_t y2=morton2 & y_mask;
    std::uint64_t z1=morton1 & z_mask; std::uint64_t z2=morton2 & z_mask;
    return x1<=x2 && y1<=y2 && z1<=z2;
}

class simple_check{
public:
    simple_check(std::uint64_t qmin, std::uint64_t qmax):query_min(qmin),query_max(qmax){}
    bool check_intersect(std::uint64_t lower_bound, std::uint64_t upper_bound) const{
        return is_less_than_or_equal(query_min,upper_bound) && is_less_than_or_equal(lower_bound,query_max);
    }
    bool check_inside(std::uint64_t morton_code) const{
        return is_less_than_or_equal(query_min,morton_code) && is_less_than_or_equal(morton_code,query_max);
    }
private:
    std::uint64_t query_min; std::uint64_t query_max;
};

class check_distance{
public:
    check_distance(const std::array<double,3>& point,double radius,int num_bits):query_point(point),radius(radius),num_bits_per_dim(num_bits){
        query_start_x=query_point[0]-radius; query_end_x=query_point[0]+radius;
        query_start_y=query_point[1]-radius; query_end_y=query_point[1]+radius;
        query_start_z=query_point[2]-radius; query_end_z=query_point[2]+radius;
        cell_shift_x_min=std::floor(query_start_x); cell_shift_x_max=std::floor(query_end_x);
        cell_shift_y_min=std::floor(query_start_y); cell_shift_y_max=std::floor(query_end_y);
        cell_shift_z_min=std::floor(query_start_z); cell_shift_z_max=std::floor(query_end_z);
    }

    bool check_point(const std::array<double,3>& p) const{
        double dx=std::abs(p[0]-query_point[0]);
        double dy=std::abs(p[1]-query_point[1]);
        double dz=std::abs(p[2]-query_point[2]);
        double dist2=dx*dx+dy*dy+dz*dz; return dist2<=radius*radius;
    }

    bool check_intersect(std::uint64_t lower_bound, std::uint64_t upper_bound) const{
        auto [xmin,ymin,zmin] = deinterleave_bits(lower_bound,num_bits_per_dim);
        auto [xmax,ymax,zmax] = deinterleave_bits(upper_bound,num_bits_per_dim);
        int level=num_bits_per_dim-1; double cell_dis=1.0/(1<<level);
        double xmin_f=morton_code_to_coordinate(xmin,num_bits_per_dim);
        double xmax_f=morton_code_to_coordinate(xmax,num_bits_per_dim)+cell_dis;
        double ymin_f=morton_code_to_coordinate(ymin,num_bits_per_dim);
        double ymax_f=morton_code_to_coordinate(ymax,num_bits_per_dim)+cell_dis;
        double zmin_f=morton_code_to_coordinate(zmin,num_bits_per_dim);
        double zmax_f=morton_code_to_coordinate(zmax,num_bits_per_dim)+cell_dis;
        if(xmin_f<=query_point[0] && xmax_f>=query_point[0] &&
           ymin_f<=query_point[1] && ymax_f>=query_point[1] &&
           zmin_f<=query_point[2] && zmax_f>=query_point[2]) return true;
        for(int sz=cell_shift_z_min; sz<=cell_shift_z_max; ++sz){
            for(int sy=cell_shift_y_min; sy<=cell_shift_y_max; ++sy){
                for(int sx=cell_shift_x_min; sx<=cell_shift_x_max; ++sx){
                    double xmin_s=xmin_f+sx; double ymin_s=ymin_f+sy; double zmin_s=zmin_f+sz;
                    double xmax_s=xmax_f+sx; double ymax_s=ymax_f+sy; double zmax_s=zmax_f+sz;
                    double dx=std::min(std::abs(xmin_s-query_point[0]), std::abs(xmax_s-query_point[0]));
                    double dy=std::min(std::abs(ymin_s-query_point[1]), std::abs(ymax_s-query_point[1]));
                    double dz=std::min(std::abs(zmin_s-query_point[2]), std::abs(zmax_s-query_point[2]));
                    double dist2=dx*dx+dy*dy+dz*dz; if(dist2<=radius*radius) return true;
                }
            }
        }
        return false;
    }

    std::vector<std::array<int,3>> check_inside(std::uint64_t morton_code) const{
        auto [xmin,ymin,zmin]=deinterleave_bits(morton_code,num_bits_per_dim);
        int level=num_bits_per_dim-1; double cell_dis=1.0/(1<<level);
        double xmin_f=morton_code_to_coordinate(xmin,num_bits_per_dim);
        double ymin_f=morton_code_to_coordinate(ymin,num_bits_per_dim);
        double zmin_f=morton_code_to_coordinate(zmin,num_bits_per_dim);
        double xmax_f=xmin_f+cell_dis; double ymax_f=ymin_f+cell_dis; double zmax_f=zmin_f+cell_dis;
        std::vector<std::array<int,3>> result;
        for(int sz=cell_shift_z_min; sz<=cell_shift_z_max; ++sz){
            for(int sy=cell_shift_y_min; sy<=cell_shift_y_max; ++sy){
                for(int sx=cell_shift_x_min; sx<=cell_shift_x_max; ++sx){
                    double xmin_s=xmin_f+sx; double ymin_s=ymin_f+sy; double zmin_s=zmin_f+sz;
                    double xmax_s=xmax_f+sx; double ymax_s=ymax_f+sy; double zmax_s=zmax_f+sz;
                    double dx=std::min(std::abs(xmin_s-query_point[0]), std::abs(xmax_s-query_point[0]));
                    double dy=std::min(std::abs(ymin_s-query_point[1]), std::abs(ymax_s-query_point[1]));
                    double dz=std::min(std::abs(zmin_s-query_point[2]), std::abs(zmax_s-query_point[2]));
                    double dist2=dx*dx+dy*dy+dz*dz; if(dist2<=radius*radius){ result.push_back({sx,sy,sz}); }
                }
            }
        }
        return result;
    }
private:
    std::array<double,3> query_point; double radius; int num_bits_per_dim;
    double query_start_x,query_end_x,query_start_y,query_end_y,query_start_z,query_end_z;
    int cell_shift_x_min,cell_shift_x_max,cell_shift_y_min,cell_shift_y_max,cell_shift_z_min,cell_shift_z_max;
};

#endif

import math 
from morton_code import *
from binary_radix_tree import *

def is_less_than_or_equal(morton1: int, morton2: int) -> bool:
    """
    Check if morton1 is less than or equal to morton2.
    """
    x_mask = 0x9249249249249249  # Mask for x bits
    y_mask = 0x2492492492492492  # Mask for y bits
    z_mask = 0x4924924924924924  # Mask for z bits
    x1 = morton1 & x_mask
    x2 = morton2 & x_mask
    y1 = morton1 & y_mask
    y2 = morton2 & y_mask
    z1 = morton1 & z_mask
    z2 = morton2 & z_mask
    return x1 <= x2 and y1 <= y2 and z1 <= z2

class simple_check:
    def __init__(self, query_min, query_max):
        self.query_min = query_min
        self.query_max = query_max
        
    def check_intersect(self, lower_bound, upper_bound):
        # Check intersection for each dimension
        # return all(query_max[i] >= lower_bound[i] and query_min[i] <= upper_bound[i] for i in range(3))
        return is_less_than_or_equal(self.query_min, upper_bound) and is_less_than_or_equal(lower_bound, self.query_max)
    
    def check_inside(self, morton_code):
        return is_less_than_or_equal(self.query_min, morton_code) and is_less_than_or_equal(morton_code, self.query_max)
    
class check_distance:
    def __init__(self, query_point, radius, num_bits_per_dim):
        self.query_point = query_point
        self.radius = radius
        self.num_bits_per_dim = num_bits_per_dim

        self.query_start_x = self.query_point[0] - self.radius
        self.query_end_x = self.query_point[0] + self.radius
        self.query_start_y = self.query_point[1] - self.radius
        self.query_end_y = self.query_point[1] + self.radius
        self.query_start_z = self.query_point[2] - self.radius
        self.query_end_z = self.query_point[2] + self.radius

        self.cell_shift_x_min = math.floor(self.query_start_x)
        self.cell_shift_x_max = math.floor(self.query_end_x)
        self.cell_shift_y_min = math.floor(self.query_start_y)
        self.cell_shift_y_max = math.floor(self.query_end_y)
        self.cell_shift_z_min = math.floor(self.query_start_z)
        self.cell_shift_z_max = math.floor(self.query_end_z)

    def check_point(self, point):
        diff_x = abs(point[0] - self.query_point[0])
        diff_y = abs(point[1] - self.query_point[1])
        diff_z = abs(point[2] - self.query_point[2])

        dist2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z
        return dist2 <= self.radius * self.radius

        
    def check_intersect(self, lower_bound, upper_bound):
        xmin, ymin, zmin = deinterleave_bits(lower_bound, self.num_bits_per_dim)
        xmax, ymax, zmax = deinterleave_bits(upper_bound, self.num_bits_per_dim)
        
        level = self.num_bits_per_dim - 1
        cell_dis = 1.0 / (1 << level)
        # print(f"cell_dis: {cell_dis}")

        xmin_f = morton_code_to_coordinate(xmin, self.num_bits_per_dim)
        xmax_f = morton_code_to_coordinate(xmax, self.num_bits_per_dim) + cell_dis
        ymin_f = morton_code_to_coordinate(ymin, self.num_bits_per_dim)
        ymax_f = morton_code_to_coordinate(ymax, self.num_bits_per_dim) + cell_dis
        zmin_f = morton_code_to_coordinate(zmin, self.num_bits_per_dim) 
        zmax_f = morton_code_to_coordinate(zmax, self.num_bits_per_dim) + cell_dis

        if (xmin_f <= self.query_point[0] and xmax_f >= self.query_point[0]
            and ymin_f <= self.query_point[1] and ymax_f >= self.query_point[1]
            and zmin_f <= self.query_point[2] and zmax_f >= self.query_point[2]):
            return True

        for shift_z in range(self.cell_shift_z_min, self.cell_shift_z_max + 1):
            for shift_y in range(self.cell_shift_y_min, self.cell_shift_y_max + 1):
                for shift_x in range(self.cell_shift_x_min, self.cell_shift_x_max + 1):
                    xmin_f_shifted = xmin_f + shift_x
                    ymin_f_shifted = ymin_f + shift_y
                    zmin_f_shifted = zmin_f + shift_z

                    xmax_f_shifted = xmax_f + shift_x
                    ymax_f_shifted = ymax_f + shift_y
                    zmax_f_shifted = zmax_f + shift_z

                    diff_x = min(abs(xmin_f_shifted - self.query_point[0]), abs(xmax_f_shifted - self.query_point[0]))
                    diff_y = min(abs(ymin_f_shifted - self.query_point[1]), abs(ymax_f_shifted - self.query_point[1]))
                    diff_z = min(abs(zmin_f_shifted - self.query_point[2]), abs(zmax_f_shifted - self.query_point[2]))

                    dist2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z
                    if dist2 <= self.radius * self.radius:
                        return True
        return False
    
    def check_inside(self, morton_code):
        xmin, ymin, zmin = deinterleave_bits(morton_code, self.num_bits_per_dim)        
        
        level = self.num_bits_per_dim - 1
        cell_dis = 1.0 / (1 << level)

        xmin_f = morton_code_to_coordinate(xmin, self.num_bits_per_dim)
        ymin_f = morton_code_to_coordinate(ymin, self.num_bits_per_dim)
        zmin_f = morton_code_to_coordinate(zmin, self.num_bits_per_dim)
        xmax_f = xmin_f + cell_dis
        ymax_f = ymin_f + cell_dis
        zmax_f = zmin_f + cell_dis
        
        result = []

        for shift_z in range(self.cell_shift_z_min, self.cell_shift_z_max + 1):
            for shift_y in range(self.cell_shift_y_min, self.cell_shift_y_max + 1):
                for shift_x in range(self.cell_shift_x_min, self.cell_shift_x_max + 1):
                    xmin_f_shifted = xmin_f + shift_x
                    ymin_f_shifted = ymin_f + shift_y
                    zmin_f_shifted = zmin_f + shift_z

                    xmax_f_shifted = xmax_f + shift_x
                    ymax_f_shifted = ymax_f + shift_y
                    zmax_f_shifted = zmax_f + shift_z

                    diff_x = min(abs(xmin_f_shifted - self.query_point[0]), abs(xmax_f_shifted - self.query_point[0]))
                    diff_y = min(abs(ymin_f_shifted - self.query_point[1]), abs(ymax_f_shifted - self.query_point[1]))
                    diff_z = min(abs(zmin_f_shifted - self.query_point[2]), abs(zmax_f_shifted - self.query_point[2]))

                    dist2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z
                    if dist2 <= self.radius * self.radius:
                        result.append([shift_x, shift_y, shift_z])
        return result
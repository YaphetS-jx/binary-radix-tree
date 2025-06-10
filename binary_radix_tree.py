import numpy as np

def count_leading_zeros(x, num_bits=30):
    return num_bits - x.bit_length() if x != 0 else num_bits

def delta(morton_codes, i, j, num_bits=30):
    n = len(morton_codes)
    if j < 0 or j >= n:
        return -1  # Sentinel for out-of-range
    return count_leading_zeros(morton_codes[i] ^ morton_codes[j], num_bits)

def determine_range(morton_codes, i, num_bits=30):
    n = len(morton_codes)
    d = np.sign(delta(morton_codes, i, i+1, num_bits) - delta(morton_codes, i, i-1, num_bits))
    delta_min = delta(morton_codes, i, i - d, num_bits)

    # Exponential search to find other end
    lmax = 2
    while delta(morton_codes, i, i + lmax * d, num_bits) > delta_min:
        lmax *= 2

    # Binary search to refine
    l = 0
    t = lmax // 2
    while t >= 1:
        if delta(morton_codes, i, i + (l + t) * d, num_bits) > delta_min:
            l += t
        t //= 2

    j = i + l * d
    return (min(i, j), max(i, j))

def find_split(morton_codes, delta_node, first, last, num_bits=30):
    split = first
    stride = last - first

    while stride > 1:
        stride = (stride + 1) // 2
        mid = split + stride
        if mid < last and delta(morton_codes, first, mid, num_bits) > delta_node:
            split = mid
    return split

def find_range_3d(prefix, num_bits_prefix, num_bits = 30):
    assert(num_bits % 3 == 0)

    # create masks 
    mask1 = ((1 << num_bits_prefix) - 1) << (num_bits - num_bits_prefix)
    mask2 = mask1 ^ ((1 << num_bits) - 1)
    
    # Create the two numbers
    lower_bound = prefix & mask1
    upper_bound = prefix | mask2

    # xmin, ymin, zmin = deinterleave_bits(min_val, num_bits//3)
    # xmax, ymax, zmax = deinterleave_bits(max_val, num_bits//3)

    # print(f"xmin: {xmin}, ymin: {ymin}, zmin: {zmin}, xmax: {xmax}, ymax: {ymax}, zmax: {zmax}")
    # return (xmin, ymin, zmin), (xmax, ymax, zmax)
    return lower_bound, upper_bound

class Node:
    def __init__(self, left, right, lower_bound, upper_bound):        
        self.left = left
        self.right = right
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

class BinaryRadixTree:
    def __init__(self, morton_codes, num_bits = 30):
        self.morton_codes = morton_codes
        self.num_bits = num_bits
        self.nodes = self.build_tree(morton_codes, num_bits)

    def build_tree(self, morton_codes, num_bits = 30):
        n = len(morton_codes)
        nodes = [None] * (n - 1)

        for i in range(n - 1):
            first, last = determine_range(morton_codes, i, num_bits)
            delta_node = delta(morton_codes, first, last, num_bits)
            split = find_split(morton_codes, delta_node, first, last, num_bits)        
            lower_bound, upper_bound = find_range_3d(morton_codes[split], delta_node, num_bits)

            left = split if split == first else split + n
            right = split + 1 if split + 1 == last else split + 1 + n
            nodes[i] = Node(
                left=left,
                right=right,
                lower_bound=lower_bound,
                upper_bound=upper_bound
            )
        return nodes

    def traverse_tree(self, check):
        results = []
        n = len(self.morton_codes)
        stack = []  # Stack to store nodes to visit
        
        # Start with root node    
        stack.append(n)
        
        while stack:
            node_idx = stack.pop()        
            
            # Handle leaf node
            if node_idx < n:
                if check.check_inside(self.morton_codes[node_idx]):
                    results.append(node_idx)
            else:
                node = self.nodes[node_idx - n]
                if (check.check_intersect(node.lower_bound, node.upper_bound)):
                    stack.append(node.left)
                    stack.append(node.right)
        return sorted(results)
    
    def traverse_tree_distance(self, check):
        results = {}
        n = len(self.morton_codes)
        stack = []  # Stack to store nodes to visit
        
        # Start with root node    
        stack.append(n)
        
        while stack:
            node_idx = stack.pop()        
            
            # Handle leaf node
            if node_idx < n:
                temp = check.check_inside(self.morton_codes[node_idx])
                if len(temp) > 0:
                    results[node_idx] = temp
                    # print(f"node_idx: {node_idx}, temp: {temp}")
                    
            else:
                node = self.nodes[node_idx - n]
                if (check.check_intersect(node.lower_bound, node.upper_bound)):
                    stack.append(node.left)
                    stack.append(node.right)
        return results

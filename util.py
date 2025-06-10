import math 

# util functions 
def prefix_sum(arr):
    prefix_sum = [0]
    for i in range(len(arr)):
        prefix_sum.append(prefix_sum[i] + arr[i])
    return prefix_sum
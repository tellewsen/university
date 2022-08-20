from numpy import array

"""
sums 1D array ARR (but IDL's total is faster and more general)
"""
ar = array([2])


def add_up(arr):
    sum_arr = 0
    for i in range(len(arr)):
        sum_arr += arr[i]
    return sum_arr


print add_up(ar)

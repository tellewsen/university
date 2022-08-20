from numpy import linspace
w = linspace(0, 3, 31)

# w[from_index:to_index_plus1:incr]
print w[:]        # w[0], w[1], ... w[20]
print w[:-2]      # w[0], w[1], ..., w[18]
print w[::5]      # w[0], w[5], w[10], w[15], w[20]
print w[2:-2:6]   # w[2], w[8], w[14]




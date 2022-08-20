This file describes the different methods this program can use

Methods:

Method 0: Print analytical solutions for expectation values for L=2,T=1 to terminal, and write computations of the same values to file
Note however that it is not really the same values, but since T= 1 it turns out that the variance of the energy is the same as the spesific heat, and this is also the case for the variance of the magnetization and the susceptibility

Method 1: Plots the number of accepted moves for each temperature in the temperature range set by the user. Here the first column is the temperature while the second is the number of accepted moves. Note that the user must keep track of how many Monte Carclo cycles and spins that were used for this computation if he or she wants to normalize their plots of this. For reference the total number of moves is the number of monte carlo cycles times the number of spins squared. (M_total = mcs*n_spin*n_spin)

Method 2: This method plots the expectation values that have been reached after the program has run ist total amount of monte carlo cycles for each temperature set by the user. In this case the columns are from left to right: temperature, expectation value of E, variance of E, expectation value of M, variance of M, and expectation value of |M|.

Method 3: This plots the expectation values computed so far for each 1000 steps through the Monte Carlo cycles for one given temperature. The columns of the output file are the same as in method 2, but the number of cycles is output at the start of each line.

Method 4: Prints a list of the energy for each Monte Carlo cycles for a given temperature. This is useful for plotting the probability of being in a given state. This output file only has one column containing the energy for each step.


Ordered vs disordered initial state:

The program lets the user choose wether he or she want an ordered or disordered initial state of the system.
this is simply control by setting the "ordered parameter to true or false"




Plotting:
There are 4 files included for plotting different things. Each file corresponds to one method, and they are named plotting1,2,3,4 where the number refers to which method was used to make the data you want to plot.
The files are python files, and you should edit the file such that the line with 'output.dat' points to the data file you want to use.

#Exercise 3.2 p. 120

def sum_1_div_k(M):  #Define function
    s=0              #Startvalue of s
    k=1              #Startvalue of k
    while k <= M:    #Checks if k has reached M
        s += 1./k    #Adding next number 
        k += 1       #Increases k by 1
    return s         #Returns the sum

print sum_1_div_k(3) #Prints the sum of 1/k from k=1 to 3



"""
Terminal>  python sum_func.py 
1.83333333333
"""


                  

initial_amount = 100
p = 5.5 # interest rate
amount = initial_amount
years = 0

while amount <= 1.5*initial_amount:
    amount += p/100.*amount
    years += 1
print "It takes %g years with %.1f%% interest rate to increase your money from %.2f to %.2f" %(years,p,initial_amount, amount)


"""
Computes how many years it takes for the initial amount to increase by 50% with interest p
"""

"""
Kjoreeksempel> python interest_rate_loop.py 
It takes 8 years with 5.5% interest rate to increase your money from 100.00 to 153.47
"""

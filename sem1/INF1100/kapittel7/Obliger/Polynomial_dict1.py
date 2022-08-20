#Exercise 7.34 p. 405

class Polynomial:

    def __init__(self, coefficients):
        self.coeff = coefficients

    def __call__(self, x):
        s = 0
        for i in self.coeff:
            s += self.coeff[i]*x**i
        return s

    def __add__(self,other):
        # Start with the longest list and add in the other
        if len(self.coeff) > len(other.coeff):
            result_coeff = dict(self.coeff) #copy!
            for i in other.coeff:
                if result_coeff[i] + other_coeff[i] != 0:
                    result_coeff[i] += other.coeff[i]
                else:
                    del result_coeff[i]
        else:
            result_coeff = dict(other.coeff) # copy!
            for i in self.coeff:
                if result_coeff[i] + self.coeff[i] != 0:
                    result_coeff[i] += self.coeff[i]
                else:
                    del result_coeff[i]

        return Polynomial(result_coeff)


#Test part(only runs if this file is run directly)
if __name__=='__main__':
    a1= Polynomial( {1:1, 100:-3} )
    a2= Polynomial( {1:-1, 20:1, 100:4} )
    a3= a1 + a2
    print a3.coeff 

"""
Terminal> python Polynomial_dict1.py
{20: 1, 100: 1}
"""

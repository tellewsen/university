#Exercise 3.9 p. 122

def area(vertices):                                             #Defines function (requires a nested list of the form [[a,b],[c,d],[e,f]] )
    x1,x2,x3=vertices[0][0],vertices[1][0],vertices[2][0]       #Makes the equation below easier to write and understand
    y1,y2,y3=vertices[0][1],vertices[1][1],vertices[2][1]       #Makes the equation below easier to write and understand

    A = 0.5*abs(x2*y3 - x3*y2 - x1*y3 + x3*y1 + x1*y2 - x2*y1)  #Calculates A (Area of a triangle)
    return A                                                    #Returns answer of A

print 'Area of triangle1 is: ', area([[0,0], [1,0], [0,1]])     #Prints the area of a triangle with vertex coords (0,0), (1,0) and (0,1). Answer should be 0.5
print 'Area of triangle2 is: ', area([[1,1], [1,4], [4,1]])     #Prints the area of a triangle with vertex coords (1,1), (1,4) and (4,1). Answer should be 4.5

"""
Terminal>  python area_triangle.py 
Area of triangle1 is:  0.5
Area of triangle2 is:  4.5
"""

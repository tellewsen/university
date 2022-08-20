def area(vertices):
    x1,x2,x3 = vertices[0][0],vertices[1][0],vertices[2][0]
    y1,y2,y3 = vertices[0][1],vertices[1][1],vertices[2][1]
    A = 0.5*abs(x2*y3 - x3*y2 - x1*y3 + x3*y1 + x1*y2 -x2*y1)
    return A


triangle1 = area([[0,0], [1,0], [0,2]])
print 'Area of triangle is %.2f' % triangle1

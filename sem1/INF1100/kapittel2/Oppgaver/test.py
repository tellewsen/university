C= []
Cd = -50
C_max= 200
while Cd <=C_max:
    C.append(Cd)
    Cd+=2.5
print C


degrees = [0, 10, 20, 40, 100]
for C in degrees:
    print "list element:", C
print "The degrees list has", len(degrees), "elements"

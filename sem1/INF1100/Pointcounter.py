"""
Input the number of points you get for each week,
and this program will tell you when you've gotten enough points.
"""



#Dictionary
weekpoints = 	{'week01':2, 'week02':3, 'week03':5, 'week04':4,
		 'week05':3, 'week06':4, 'week07':5, 'week08':3,
		 'week09':2, 'week10':4, 'week11':4, 'week12':0}
totalpoints = 0

#for-loop
for i in weekpoints:
    totalpoints += weekpoints[i]

#if-test
if totalpoints >= 35:
    print "You now have %g points, that's enough!" %(totalpoints)
    print "You're AWESOME!!!!!!!!!!!!!!!!!!!!!"
else:
    print "You now have %g points, you need 35!" %(totalpoints)
    print "Get to WORK!!!!"

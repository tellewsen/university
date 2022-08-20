import calendar,sys 
a = calendar.weekday(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])) #yyyy,mm,dd
print calendar.day_name[a] #prints the name of the day you input on the cml


"""
Terminal>  python weekday.py 1988 1 4
Monday

Terminal> 1988 1 5
Tuesday

Terminal> python weekday.py 1989 1 4
Wednesday

Terminal> weekday.py 1989 7 30
Sunday
"""

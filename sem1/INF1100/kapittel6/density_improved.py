#Exercise 6.10 p.333

def read_densities1(filename):		#Part 1 of exercise
    infile = open(filename, 'r')	#Opens file
    densities = {}			#Makes dictionary densities
    for line in infile:			
        words = line.split()		#Splits words
        density = float(words[-1])	#Converts last part into float and calls it density
	substance = ' '.join(words[:-1])#Puts the rest of the line in substance
        densities[substance] = density	#Stores density inside substance inside densities
    infile.close()			#Close file
    return densities			


def read_densities2(filename):		#Part 2 of exercise
    infile = open(filename, 'r')	#Opens file
    densities = {}			#Makes dictionary densities
    for line in infile:
        density = float(line[12:])	#Makes float of character 12 to end of line and store in density
        substance = line[:12].rstrip()	#From start to character 12, then strip spaces at end.
        densities[substance] = density	#Stores density inside substance inside densities
    infile.close()			#Close file
    return densities



#Tests that you actually get the same thing out of both functions
if read_densities1('densities.dat') == read_densities2('densities.dat'):
    print 'both functions print the same thing'

"""
Terminal>  python density_improved.py 
both functions print the same thing
"""

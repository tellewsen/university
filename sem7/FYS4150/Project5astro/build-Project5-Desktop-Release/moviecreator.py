from scitools.easyviz import *
import glob,os
movie('tmp*.png',encoder='mencoder',fps = 24)

for filename in glob.glob('tmp*.png'):
	os.remove(filename)

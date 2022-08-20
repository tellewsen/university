;Lab 2 for AST2210
;Author: Andreas Ellewsen

;Exercise7
bias1 = read_bmp('bias1.bmp')
bias2 = read_bmp('bias2.bmp')
biassum = bias1 + bias2
biasdiff = bias1 - bias2
biassquare1 = biassum[376-149:376+150,240-149:240+150]
biassquare2 = biasdiff[376-149:376+150,240-149:240+150]
meanbiassquare1 = mean(biassquare1)
stddevbiassquare2 = stddev(biassquare2)

print,'Mean of central region sum: ' + strtrim(meanbiassquare1,2)
print,'Std dev of central region diff: ' + strtrim(stddevbiassquare2,2)

;Exercise9
flat3 = read_bmp('flat3.bmp')
flat4 = read_bmp('flat4.bmp')
flatsum  = flat3 + flat4
flatdiff = flat3 - flat4
flatsquare1 = flatsum[376-149:376+150,240-149:240+150]
flatsquare2 = flatdiff[376-149:376+150,240-149:240+150]
meanflatsquare1 = mean(flatsquare1)
stddevflatsquare2 = stddev(flatsquare2)

print,'Mean of central region sum: ' + strtrim(meanflatsquare1,2)
print,'Std dev of central region diff: ' + strtrim(stddevflatsquare2,2)

;Exercise10
g = (meanflatsquare1 - meanbiassquare1) / (stddevflatsquare2^2 - stddevbiassquare2^2)
print,'Conversion factor g = ',strtrim(g,2)

;Exercise 11 of lab2 for AST2210
readoutnoise = stddevbiassquare2*g
print,'Readout noise for bias frames = ',strtrim(readoutnoise,2)

;Exercise12 of lab2 for AST2210
flat5  = read_bmp('flat5.bmp')
flat6  = read_bmp('flat6.bmp')
flat7  = read_bmp('flat7.bmp')
flat8  = read_bmp('flat8.bmp')
flat9  = read_bmp('flat9.bmp')
flat10 = read_bmp('flat10.bmp')
flat11 = read_bmp('flat11.bmp')
flat12 = read_bmp('flat12.bmp')
flat13 = read_bmp('flat13.bmp')
flat14 = read_bmp('flat14.bmp')
flat15 = read_bmp('flat15.bmp')
flat16 = read_bmp('flat16.bmp')

stddevflatsquarefinal = make_array(1,7,/float,value = 0)
flatdiff = flat3 - flat4
flatsquarefinal = (flatdiff[376-149:376+150,240-149:240+150])/sqrt(2)
stddevflatsquarefinal[0] = stddev(flatsquarefinal)

flatdiff = (flat3+flat5) - (flat4+flat6)
flatsquarefinal = (flatdiff[376-149:376+150,240-149:240+150])/sqrt(4)
stddevflatsquarefinal[1] = stddev(flatsquarefinal)

flatdiff = (flat3+flat5+flat7) - (flat4+flat6+flat8)
flatsquarefinal = (flatdiff[376-149:376+150,240-149:240+150])/sqrt(6)
stddevflatsquarefinal[2] = stddev(flatsquarefinal)

flatdiff = (flat3+flat5+flat7+flat9) - (flat4+flat6+flat8+flat10)
flatsquarefinal = flatdiff[376-149:376+150,240-149:240+150]/sqrt(8)
stddevflatsquarefinal[3] = stddev(flatsquarefinal)

flatdiff = (flat3+flat5+flat7+flat9+flat11) - (flat4+flat6+flat8+flat10+flat12)
flatsquarefinal = flatdiff[376-149:376+150,240-149:240+150]/sqrt(10)
stddevflatsquarefinal[4] = stddev(flatsquarefinal)

flatdiff = (flat3+flat5+flat7+flat9+flat11+flat13) - (flat4+flat6+flat8+flat10+flat12+flat14)
flatsquarefinal = flatdiff[376-149:376+150,240-149:240+150]/sqrt(12)
stddevflatsquarefinal[5] = stddev(flatsquarefinal)

flatdiff = (flat3+flat5+flat7+flat9+flat11+flat13+flat15) - (flat4+flat6+flat8+flat10+flat12+flat14+flat16)
flatsquarefinal = flatdiff[376-149:376+150,240-149:240+150]/sqrt(14)
stddevflatsquarefinal[6] = stddev(flatsquarefinal)

window,0
plot,[0,1,2,3,4,5,6],stddevflatsquarefinal, title='Noise in flats',xtitle= 'Number of flats used *2 +2',ytitle = 'Noise[relevant units]'

write_png,'Noiseinflats.png',tvrd()
;As expected we see a decrease in the noise for every increase in
;number of flats used. If any programmers come across this code, I
;apologize for the horrible piece at the end here.

;Exercise13

Rawimage = read_bmp('CCD single slit.bmp')
darkflat1 = read_bmp('darkflat1.bmp')
darkflat2 = read_bmp('darkflat2.bmp')
darkflat3 = read_bmp('darkflat3.bmp')
darkflat4 = read_bmp('darkflat4.bmp')
darkflat5 = read_bmp('darkflat5.bmp')
dark1 = read_bmp('dark1.bmp')
dark2 = read_bmp('dark2.bmp')
dark3 = read_bmp('dark3.bmp')
dark4 = read_bmp('dark4.bmp')
dark5 = read_bmp('dark5.bmp')

Df = (darkflat1+darkflat2+darkflat3+darkflat4+darkflat5)/5.0 ;Master darkframe
Ddiff = (dark1+dark2+dark3+dark4+dark5)/5.0
F = (flat3+flat4+flat5+flat6+flat7+flat8+flat9+flat10+flat11+flat12+flat13+flat14+flat15+flat16)/14.0 ;Master flatfield
m = avg(F-Df)
Corrected = (Rawimage-Ddiff)*m /(F-Ddiff)

window,1
plot_image,Rawimage
write_png,'Raw image.png',tvrd()
window,2
plot_image,Corrected
write_png,'Corrected.png',tvrd()
end

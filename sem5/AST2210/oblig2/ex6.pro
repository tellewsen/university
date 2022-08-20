bias1   = read_bmp('bias1.bmp')
darkmax = read_bmp('darkmax.bmp')

;Print information for each image to terminal
print,'Info for bias1'
print,'Maximum:  '       + strtrim(max(fix(bias1)),2)
print,'Minimum:  '       + strtrim(min(fix(bias1)),2)
print,'Average:  '       + strtrim(avg(fix(bias1)),2)

print,'Info for darkmax'
print,'Maximum:  '       + strtrim(max(fix(darkmax)),2)
print,'Minimum:  '       + strtrim(min(fix(darkmax)),2)
print,'Average:  '       + strtrim(avg(fix(darkmax)),2)

window,0
plot,histogram(bias1),xrange = [0,20], title='Bias frame 1', xtitle = 'Pixel value', ytitle = 'Number of pixels with pixel value'
write_png,'histogram1.png',tvrd()
window,1
plot,histogram(darkmax),xrange = [0,20],title='Dark frame max',xtitle = 'Pixel value', ytitle = 'Number of pixels with pixel value'
write_png,'histogram2.png',tvrd()

maxbias1 = where(bias1 eq max(bias1))
maxplacementbias1 = array_indices(bias1,maxbias1)
print,maxplacementbias1

maxdarkmax = where(darkmax eq max(darkmax))
maxplacementdarkmax = array_indices(darkmax,maxdarkmax)
print,maxplacementdarkmax

end

;To begin with we notice that all values are higher for the maximum
;exposure dark frame. However, the difference is not by much, except
;for the maximum. As we see in the table, the max for the dark frame
;is 94, while the max for the bias frame is only 16.

;Locating the coordinates of the maximums reveals that they are not
;the same. This is not as expected.
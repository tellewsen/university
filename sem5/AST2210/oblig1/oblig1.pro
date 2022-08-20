;Oblig1 AST2210
;This assignment goes through how to read images from the
;hinodedata. 

;Exercise 3

;Read files into IDL
o=obj_new('osdc')
o->show,'INSTRUME,DATE_OBS,SX__PROG_VER,THUMBS'
o->page,1
o->limit,50 ;; Use o->limit,0 to get all!
o->order,'DATE_OBS',/DESCENDING
o->condition,'INSTRUME: ,SOT/WB'
o->condition,'EPOCH_START: 01-Apr-2007'
o->condition,'EPOCH_END: 30-Apr-2007'
o->condition,'SS_L1LEAD: y'
o->condition,'SX__PROG_VER: 52'
o->condition,'SQ_OBSTITLE_PROG_VER: FG20070429_182200.6'
o->search,out
paths = o->paths(out,/fetch) ;; /unzip has been deprecated
paths=paths[1:*]
read_sot,paths,index,data

;Plot the images, set title of each of them, then save
for i=0,5 do begin
image = data[*,*,i]
window,i
savep,'file'+strtrim(i,2)
plot_image,image,title=index[i].wave
savep

;Print information for each image to terminal
print,'Info for image '  + strtrim(index[i].wave)
print,'Maximum:  '       + strtrim(max(image),2)
print,'Minimum:  '       + strtrim(min(image),2)
print,'Mean:     '       + strtrim(mean(image),2)
print,'Std dev:  '       + strtrim(stddev(image),2)
print,' '

endfor

end


;for i=0,5 do begin
;
;;Exercise 4
;;Fix irregularities 
;fg_prep,index[i],data[*,*,i],index2,data2,dark_image=dark_image,flat_image=flat_image
;window,6 + 2*i
;plot_image,dark_image
;window,6 + 2*i+1
;plot_image,flat_image
;endfor



;0: flat_image has "currents" going around, there's a vertical
;line in the middle(the line between the two CCDs). dark_image has one light
;part to the left and a darker(still pretty light) on the right.

;1: dark_image has faint, dark horizontal lines going through. flat_image
;now has one vertial line in the middle of each side.

;2: flat_image has diagonal lines going from top right to bottom
;left. The previously mentioned lines are still present. dark_image
;has the horizontal lines already mentioned as well, but they are less
;faint now.

;3: flat_image now has diagonal lines going from top left to bottom
;right. All previous lines still present. dark_image still has
;horizontal lines. Not sure if they are more visible or not.

;4: Same as before, but more visible for both images. Noticing a
;pattern emerging in flat_image. 

;5: Same thing again. Noticing small white and black dots appearing
;in flat_image. dark_image still has the horizontal lines, some more
;visible than others. 

;I'm guessing that the pattern in flat_image shows the 


;Exercise 5
;fg_prep,index,data,index2,data2

;for i = 0,5 do begin
;window,5+i
;plot_image,image,title=index[i].wave
;print,'Contrast in percent:  '       + strtrim(stddev(image)/mean(image)*100,2)
;endfor


;In the unaltered picture (red cont 6684) we can see the distortion
;from what we assume is dark current, and the difference in brightness
;between the two halves. Both of these errors are fixed in the altered
;image, and the overall brightness has been lowered. However the line
;between the two halves is still visible.

;When listing contrast in percent for the altered images we get 11.2689 for
;every image

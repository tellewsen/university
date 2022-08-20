;Oblig1 AST2210
;Exercise 6

;Read files into IDL
o=obj_new('osdc')
o->show,'FILE,INSTRUME,DATE_OBS,CEN_RADIUS,S__WAVE,THUMBS'
o->page,1
o->limit,0 ;; Use o->limit,0 to get all!
o->order,'DATE_OBS',/ASCENDING
o->condition,'INSTRUME: ,SOT/WB'
o->condition,'CEN_RADIUS: >1000'
o->condition,'SS__L1LEAD: y'
o->condition,'Gx: NONE'
o->condition,'S__WAVE_ID_: ,7'
o->condition,'xSQ_IUMODE1: FG20130930_193029.5'
o->search,out
paths =  o->paths(out,/fetch) ;; /unzip has been deprecated
paths=paths[1:*]
;read_sot,paths,index,data

;fg_prep,index,data,index2,data2

mfg_mkcube,paths,'cah'
restore,'cah.idl',/verbose
save,hdr,file='cah.idl'
crispex,'cah.icube', /single

image=lp_get('cah.icube',20) ; read image number 20 from the timeseries
end

;We do as the assignment tells us to do, choosing CEN_RADIUS, SOT
;search criterium S__WAVE with Ca II H Line. We use the example code
;given in the assignment and start the program. After about 15 minutes
;the program is done and the crispex window opens. It turns out
;I've made a rather poor choice when selecting which series of
;images to use. We do see very clearly that the limb moves throughout
;the movie.

;Exercise7:
lp_header,'cah.icube',nt=nt ; determine the number of images
dxy=fltarr(2,nt) ; make a dxy-array filled with zeroes
dx = (809-786)/100.
dy = 0
dxy[0,*] = dx ; average displacement per timestep in x direction
dxy[1,*] = dy ; average displacement per timestep in y direction
mfg_shift,'cah.icube',dxy ; compensate for the drift

end

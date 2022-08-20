;Load files into IDL
f   = file_search('ALvi1600*.fits') 
fn  = long(stregex(f,'[0-9]{6}',/ext))		
fb  = f[where(fn - 160000 lt 10)]                            ;bias
ffO = f[where((fn - 160000 ge 40) and (fn - 160000 le 43) )] ;flatfields OIII
ffH = f[where((fn - 160000 ge 44) and (fn - 160000 le 47) )] ;flatfields Ha


;Mean of bias and flatfields
for i =  0,n_elements(fb)-1 do begin
	im = readfits(fb[i],exten_no=1)
	if i eq 0 then bias = im else bias = [ [[bias]] , [[im]] ] 
     endfor

for i =  0,n_elements(ffO)-1 do begin
	im = readfits(ffO[i],exten_no=1)
	if i eq 0 then flatsO = im else flatsO = [ [[flatsO]] , [[im]] ] 
     endfor

for i =  0,n_elements(ffH)-1 do begin
	im = readfits(ffH[i],exten_no=1)
	if i eq 0 then flatsH = im else flatsH = [ [[flatsH]] , [[im]] ] 
endfor

meanffH  = float(total(flatsH,3))/n_elements(ffH)
meanffO  = float(total(flatsO,3))/n_elements(ffO)
meanbias = float(total(bias,3))/n_elements(bias)


;Read Data images

imO = f[where((fn - 160000 ge 78) and (fn - 160000 le 79) )]
imH = f[where((fn - 160000 ge 80) and (fn - 160000 le 81) )]

for i =  0,n_elements(imO)-1 do begin
	im = readfits(imO[i],exten_no=1)
	if i eq 0 then imageO = im else imageO = [ [[imageO]] , [[im]] ] 
endfor

for i =  0,n_elements(imH)-1 do begin
	im = readfits(imH[i],exten_no=1)
	if i eq 0 then imageH = im else imageH = [ [[imageH]] , [[im]] ] 
endfor

;Correct images

for i = 0,n_elements(imO)-1 do begin
       corrO = float(imageO[*,*,i]-meanbias)/(meanffO-meanbias)
       if i eq 0 then corrimO = corrO else corrimO = [ [[corrimO]] , [[corrO]] ] 
endfor

for i = 0,n_elements(imH)-1 do begin
      corrH = float(imageH[*,*,i]-meanbias)/(meanffH-meanbias)
      if i eq 0 then corrimH = corrH else corrimH = [ [[corrimH]] , [[corrH]] ] 
endfor

;plot images
;window,0
;plot_image,corrimO[*,*,0] < 0.075
;plot,corrimO[*,1000,0] < 0.075
;window,1
;plot_image,corrimO[*,*,1] < 0.075
;plot,corrimO[*,1000,1] < 0.075
;window,2
;plot_image,corrimH[*,*,0] < 0.075
;window,3
;plot_image,corrimH[*,*,1] < 0.075

;Plot images on top of each other
bestimO = corrimO[*,*,0] < 0.075
bestimH = corrimH[*,*,0] < 0.075
window,0,xs=2148/3,ys=2102/3

;Shift images to overlap better
;bestimH = shift(bestimH,[-100,-100])
;bestimH = 
tvscale,congrid(bestimH,2148/3,2102/3),channel=2
tvscale,congrid(bestimO,2148/3,2102/3),channel=3

;Get the files for the spectrums

imS = f[where( fn - 160000 eq 91)]                           ;NGC spectrum
Wss = f[where((fn - 160000 ge 76) and (fn - 160000 le 77) )] ;Wolf 1346 spectrum
Sff = f[where((fn - 160000 ge 87) and (fn - 160000 le 88) )] ;Spectroscopic flat fields
wlc = f[where((fn - 160000 ge 89) and (fn - 160000 le 90) )] ;Wavelength calibration

;Get the spectrums
for i =  0,n_elements(Wss)-1 do begin
	temp = readfits(Wss[i],exten_no=1)
	if i eq 0 then  Wolfspec = temp else Wolfspec = [ [[Wolfspec]] , [[temp]] ] 
endfor
for i =  0,n_elements(Sff)-1 do begin
	temp = readfits(Sff[i],exten_no=1)
	if i eq 0 then  Specflat = temp else Specflat = [ [[Specflat]] , [[temp]] ] 
endfor
for i =  0,n_elements(wlc)-1 do begin
	temp = readfits(wlc[i],exten_no=1)
	if i eq 0 then  wavecali = temp else wavecali = [ [[wavecali]] , [[temp]] ]
endfor

NGCspec = readfits(imS,exten_no=1)

;Construct means & correct
meanspec  = float(total(Specflat,3))/n_elements(Sff)
meanwlc   = float(total(Wavecali,3))/n_elements(wlc)
meanWolf  = float(total(Wolfspec,3))/n_elements(Wss)


;Correct the Wolf spectra
for i = 0,n_elements(Wss)-1 do begin
       temp = float(Wolfspec[*,*,i])/(meanspec)*(meanwlc)
       if i eq 0 then corrWolf = temp else corrWolf = [ [[corrWolf]] , [[temp]] ]
endfor

;Correct the NGC spectrum
for i = 0,n_elements(imS)-1 do begin
       temp = (float(NGCspec)/(meanspec))*(meanwlc)
       if i eq 0 then corrNGC = temp else corrNGC = [ [[corrNGC]] , [[temp]] ] 
endfor



;print stuff to terminal

help,Wolfspec
help,Specflat
help,wavecali
help,NGCspec

help, meanspec
help, meanwlc
help, meanWolf
help, corrWolf
help, corrNGC

;plot NGC, Wolf1 and 2
;window, 0
;plot,(total(corrNGC, 2))[600:1500], /yn

;window, 1
;plot,total(corrWolf[*,*,0], 2), /yn

;window, 3
;plot,total(corrWolf[*,*,1], 2), /yn


;window,1,xs=1600/1.15,ys=900/1.15
;plot,total(Wolfspec,1),/yn

;NGCspec = readfits(imS,h,exten_no=1)
;help,NGCspec
;window,1,xs=1600/1.15,ys=900/1.15
;plot,NGCspec

end

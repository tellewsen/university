function addup,arr
;+
; sums 1D array ARR (but IDL's total is faster and more general)
;-
arraysize=SIZE(arr)
if (arraysize[0] ne 1) then print,'addup input is not a 1D array'
sumarr=0
for i=0,arraysize[1]-1 do sumarr=sumarr+arr[i]
return,sumarr
end  

    

pro vsb, nb=nb, radix=radix, dir=dir
if not(keyword_set(nb)) then nb=0
if not(keyword_set(radix)) then radix="Density"
if not(keyword_set(dir)) then dir="/data/frederic/fargo3d/shearingboxmriu/"
number = strcompress('000000'+string(nb),/remove_all)
strnb = strmid(number,5,/reverse)
filename = dir+radix+strnb+".dat"

m=dblarr(150,150,40)

openr,1,filename
readu,1,m
close,1

erase

tvscl,m(*,*,20),180,180

;we now evaluate the amount of shift of the inner and outer boxes
; and we display them adequately shifted.

date = nb*0.314159265359

shift = 2.*0.75*4.*date
shift = shift/6.
shift = shift-long(shift)
shift = shift*150.
tvscl,m(*,*,20),180+shift,30
tvscl,m(*,*,20),30+shift,30
tvscl,m(*,*,20),180-shift,330
tvscl,m(*,*,20),330-shift,330
print,minmax(m)
return
end


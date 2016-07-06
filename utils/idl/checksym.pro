pro checksym, nb=nb, radix=radix, sign=sign, ymean=ymean, zmean=zmean
window,0
window,2
if not(keyword_set(radix)) then radix="Density"
if not(keyword_set(nb)) then nb=0
if not(keyword_set(sign)) then sign=1
filename=strcompress('000000'+string(nb),/remove_all)
filename=strmid(filename,strlen(filename)-6,6)
filename=radix+filename+".dat"
print,filename
m=dblarr(256,256)
openr,1,filename
readu,1,m
close,1
;m=shift(m,0,30)
if (keyword_set(ymean)) then m=.5*(m+shift(m,-1,0))
if (keyword_set(zmean)) then m=.5*(m+shift(m,0,-1))
wset,0
tvscl,m
d=m-sign*reverse(reverse(m,1),2)
tvscl,d,260,0
print,minmax(m),minmax(d)
tv,bytscl(d,min=-1d-8,max=1d-8),520,0
wset,2
shade_surf,d,charsize=2.5,zr=[-1,1]*1d-11
;xsurface,d
return
end

	


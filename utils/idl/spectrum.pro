pro spectrum, dir=dir, radix=radix, nb=nb, m=m, spec=spec,help=help
;+
; NAME:
;       SPECTRUM
;
; PURPOSE:
;
;       Extract spectra from data cubes provided
;       by FARGO3D
;
; CATEGORY:
;       Data analysis routine
;
; CALLING SEQUENCE:
;
;       SPECTRUM [,DIR=dir] [,NB=nb] [,SPEC=spec] [,M=m]
;       [,RADIX=radix] [,/HELP]
;
; DIR=dir: the directory where to look for the files. Defaults to './'
;
; NB=nb: the output number of the data that must be used. Defaults to
; 0
;
; RADIX=radix: a string that indicates the field of which we evaluate
; the spectrum. Defaults to "Density".
;
; M=m: upon return "m" is a 1D array containing the values of the
; azimuthal numbers. It is simply a linear ramp up to the value of the
; Nyquist wavenumber (Nx/2).
;
; SPEC=spec: upon return contains the power spectrum of the desired
; field (see RADIX), averaged in z and in radius.
;
; /HELP: prints this help.
;
;-
if keyword_set(help) then begin
    doc_library, 'spectrum'
    return
 endif

if not(keyword_set(nb)) then nb=0
if not(keyword_set(dir)) then dir="./"
if not(keyword_set(radix)) then radix="Density"


; Firstly we check the dimension of the simulation.
readcol,dir+"/dims.dat",f,f,f,f,f,f,nx,ny,nz,/silent
nx = nx[0]
ny = ny[0]
nz = nz[0]
; next we check the number of lines in domain_y:
com = "wc "+dir+"/domain_y.dat"
spawn,com,res
nyg=long(res)
ngh = (nyg-1-ny)/2 ; we infer the number of ghosts

readcol,dir+"/domain_y.dat",rad,/silent
rad = rad(ngh:ngh+ny)
rmed = .5*(rad+rad(1:*))
dr = rad(1:*)-rad

readcol,dir+"/domain_z.dat",z,/silent
z = z(ngh:ngh+nz)
zmed = .5*(z+z(1:*))
dz = z(1:*)-z

readcol, dir+"/domain_x.dat",x,/silent
Delta_phi = max(x)-min(x)

; we now dimension the arrays

r   = dblarr(nx,ny,nz)

; and we read them

strnum = strcompress('000000'+string(nb),/remove_all)
strnum = strmid (strnum, 5, /reverse)

openr,1,dir+"/"+radix+strnum+".dat"
readu,1,r
close,1

ms = fft(r, dimension=1)
ms = ms*conj(ms)
spec = total(total(ms,2),2)
spec = spec(0:nx/2-1)
m = indgen(nx/2-1)
return
end


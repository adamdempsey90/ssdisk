pro alpha, dir=dir, nb=nb, alpha=alpha, radius=radius, help=help, reynolds=reynolds, maxwell=maxwell
;+
; NAME:
;       ALPHA
;
; PURPOSE:
;
;       Extract the alpha value from data cubes provided
;       by FARGO3D
;
; CATEGORY:
;       Data analysis routine
;
; CALLING SEQUENCE:
;
;       ALPHA [,DIR=dir] [,NB=nb] [,ALPHA=alpha] [,RADIUS=radius]
;       [,/HELP] [,/MAXWELL] [,/REYNOLDS]
;
; DIR=dir: the directory where to look for the files. Defaults to './'
;
; NB=nb: the output number of the data that must be used. Defaults to
; 0
;
; /HELP: prints this help.
;
; ALPHA=alpha: upon return "alpha" is a 1D array that has one zone
; LESS than the original data cube (due to radial averaging issues),
; and which contains the estimate of the alpha coefficient
;
; RADIUS=radius: upon return "radius" is a 1D array that has one zone
; LESS than the original data cube (due to radial averaging issues),
; and which contains the radius of the zone centers.
;
; /MAXWELL: if this keyword is set, returns the alpha value arising
; from the Maxwell stress tensor only
;
; /REYNOLDS: if this keyword is set, returns the alpha value arising
; from the Reynolds stress tensor only
;
;-
if keyword_set(help) then begin
    doc_library, 'alpha'
    return
 endif

if not(keyword_set(nb)) then nb=0
if not(keyword_set(dir)) then dir="./"

MU0 = 1D0


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

dsurf = .5*Delta_phi/nx*dz(0); note: this routine should be modified for uneven mesh spacings

; we now dimension the arrays

vr   = dblarr(nx,ny,nz)
vphi = dblarr(nx,ny,nz)
br   = dblarr(nx,ny,nz)
bphi = dblarr(nx,ny,nz)
rho  = dblarr(nx,ny,nz)
cs   = dblarr(nx,ny,nz)

; and we read them

strnum = strcompress('000000'+string(nb),/remove_all)
strnum = strmid (strnum, 5, /reverse)

openr,1,dir+"/Bx"+strnum+".dat"
readu,1,bphi
close,1

openr,1,dir+"/By"+strnum+".dat"
readu,1,br
close,1

openr,1,dir+"/Vx"+strnum+".dat"
readu,1,vphi
close,1

openr,1,dir+"/Vy"+strnum+".dat"
readu,1,vr
close,1

openr,1,dir+"/Density"+strnum+".dat"
readu,1,rho
close,1

openr,1,dir+"/Energy"+strnum+".dat"
readu,1,cs
close,1

bymed = .5*(br+br(*,1:*,*))
help,bymed
bxmed = .5*(bphi+shift(bphi,1,0,0))
bxmed = bxmed(*,0:ny-2,*)

Tmaxwell = total(total(bxmed*bymed/MU0,3),1)*dsurf
vphiavg = total(total(vphi,3),1)/nz/nx
vravg = total(total(vr,3),1)/nz/nx

dvphi = vphi
dvrad = vr

vphiavgarr = vphiavg##replicate(1,nx,1)
vravgarr   = vravg##replicate(1,nx,1)

for k=0,nz-1 do begin
   dvphi(*,*,k) = vphi(*,*,k)-vphiavgarr
   dvrad(*,*,k) = vr(*,*,k)-vravgarr
endfor

dvphimed = 5D-1*(dvphi+shift(dvphi,1,0,0))
dvradmed = 5D-1*(dvrad+dvrad(*,1:*,*))
dvphimed = dvphimed(*,0:ny-2,*)
rho = rho(*,0:ny-2,*)

Treynolds = total(total(rho*dvradmed*dvphimed,3),1)*dsurf

cs2mean = total(total(cs*cs,3),1)/nx/nz
sigma   = total(total(rho,3),1)*dsurf/Delta_Phi

alpha = (treynolds-tmaxwell)/(Delta_Phi*sigma*cs2mean)

radius=rmed(0:ny-2)
return
end


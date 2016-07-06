pro mag_tilt_angle, nb=nb, dir=dir, angle=angle, pmag=pmag, help=help
;+
; NAME:
;       MAG_TILT_ANGLE
;
; PURPOSE:
;
;       Extract the magnetic tilt angle from data cubes
;       provided by FARGO3D
;
; CATEGORY:
;       Data analysis routine
;
; CALLING SEQUENCE:
;
;       MAG_TILT_ANGLE [,DIR=dir] [,NB=nb] [,ANGLE=angle] [,PMAG=pmag]
;       [,/HELP]
;
; DIR=dir: the directory where to look for the files. Defaults to './'
;
; NB=nb: the output number of the data that must be used. Defaults to
; 0
;
; ANGLE=angle: upon return a 3D array that contains the value of the
; tilt angle as defined by Guan et al. (2009) ApJ, 694, 1010. This
; array is one zone shorter in radius and in z (or colatitude)
;
; PMAG = pmag: upon return a 3D array that contains the 
; magnetic pressure (one zone shorter in radius and in z
; or colatitude, as ANGLE).
;
; /HELP: prints this help.
;
;-
if keyword_set(help) then begin
    doc_library, 'mag_tilt_angle'
    return
 endif

if not(keyword_set(nb)) then nb=0
if not(keyword_set(dir)) then dir="./"


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

bx   = dblarr(nx,ny,nz)
by   = dblarr(nx,ny,nz)
bz   = dblarr(nx,ny,nz)

; and we read them

strnum = strcompress(string(nb),/remove_all)
;strnum = strmid (strnum, 5, /reverse)

openr,1,dir+"/bx"+strnum+".dat"
readu,1,bx
close,1

openr,1,dir+"/by"+strnum+".dat"
readu,1,by
close,1

openr,1,dir+"/bz"+strnum+".dat"
readu,1,bz
close,1

;Periodic extension of the arrays
bbx = dblarr(nx+1,ny,nz)
bbx(0:nx-1,*,*) = bx
bbx(nx,*,*) = bx(0,*,*)

;zone centered averages:
bx = .5*(bbx(0:nx-1,*,*)+bbx(1:nx,*,*))
by = .5*(by(*,0:ny-2,*)+by(*,1:ny-1,*))
bz = .5*(bz(*,*,0:nz-2)+by(*,*,1:nz-1))

; truncate arrays to have same size
bx = bx(*,0:ny-2,0:nz-2)
by = by(*,0:ny-2,0:nz-2)
bz = bz(*,0:ny-2,0:nz-2)

help,bx,by,bz

pmag = (bx*bx+by*by+bz*bz)/2.

angle = asin(abs(bx*by)/pmag)*.5
return
end



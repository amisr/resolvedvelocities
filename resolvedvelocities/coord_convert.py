# coord_convert.py

# Series of routines to perform geographic coordinate conversions (NO GEOMAGNETIC)
# Equations/methods are taken from the following:
#   Laundal, K. M. and Richmond, A. D. (2017). Magnetic Coordinate Systems. 
#       Space Sci Rev, 206:27-59. doi: 10.1007/s11214-016-0275-y
#   Zhu, J. (1994). Conversion of Earth-centered Earth-fixed coordinates to geodetic coordinates. 
#       IEEE Trans Aerosp Electron Syst, 30(3): 957-961. doi: 10.1109/7.303772

# Throughout, the following naming and unit conventions are used:
#   cartesian: x (m), y (m), z (m)
#   spherical: r (m), t (rad), p (rad)
#   geocentric: gclat (deg), gclon (deg), gcalt (km)
#   geodetic: gdlat (deg), gdlon (deg), gdalt (km)
#   * the center of the earth is the orgin
#   * gcalt is from the center of the earth while gdalt is from the surface of the earth
#   * the latitude coordinate in geocentric/geodetic is the latitude (measured from equator)
#   * the theta coordinate in spherical is the colatitude (measured from noth pole)
#   * apologies for mixing m/km and rad/deg - this choice was made because most utilitities that have
#       geodetic/geocentric input seem to expect km/deg units while most utilities that have 
#       cartesian/spherical input seem to expect m/rad
#
# author - Leslie Lamarche
# updated - 2019-03-26

import numpy as np

# parameters defining the ellipsoid earth (from WGS84)
Req = 6378137.
f = 1/298.257223563
e2 = 2*f-f**2

def cartesian_to_spherical(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    t = np.arccos(z/r)
    p = np.arctan2(y,x)
    return r, t, p

def spherical_to_cartesian(r,t,p):
    x = r*np.sin(t)*np.cos(p)
    y = r*np.sin(t)*np.sin(p)
    z = r*np.cos(t)
    return x, y, z

def geocentric_to_spherical(gclat,gclon,gcalt):
    r = gcalt*1000.
    t = (90.-gclat)*np.pi/180.
    p = gclon*np.pi/180.
    return r, t, p

def spherical_to_geocentric(r,t,p):
    gclat = 90.-t*180./np.pi
    gclon = p*180./np.pi
    gcalt = r/1000.
    return gclat, gclon, gcalt

def cartesian_to_geocentric(x,y,z):
    r, t, p = cartesian_to_spherical(x, y, z)
    gclat, gclon, gcalt = spherical_to_geocentric(r,t,p)
    return gclat, gclon, gcalt

def geocentric_to_cartesian(gclat,gclon,gcalt):
    r, t, p = geocentric_to_spherical(gclat, gclon, gcalt)
    x, y, z = spherical_to_cartesian(r,t,p)
    return x, y, z

def geodetic_to_cartesian(gdlat,gdlon,gdalt):
    lam_gd = gdlat*np.pi/180.
    h = gdalt*1000.
    phi = gdlon*np.pi/180.
    rho = Req/np.sqrt(1-e2*np.sin(lam_gd)**2)
    x = (rho+h)*np.cos(lam_gd)*np.cos(phi)
    y = (rho+h)*np.cos(lam_gd)*np.sin(phi)
    z = (rho+h-e2*rho)*np.sin(lam_gd)
    return x, y, z

def cartesian_to_geodetic(x,y,z):
    # Heikkinen method
    a = Req
    b = a*np.sqrt(1-e2)
    r = np.sqrt(x**2+y**2)
    a2 = a**2
    b2 = b**2
    r2 = r**2
    z2 = z**2
    ep2 = (a2-b2)/b2
    F = 54*b2*z2
    G = r2+(1-e2)*z2-e2*(a2-b2)
    c = e2**2*F*r2/G**3
    s = np.cbrt(1+c+np.sqrt(c**2+2*c))
    P = F/(3*(s+1/s+1)**2*G**2)
    Q = np.sqrt(1+2*e2**2*P)
    r0 = -P*e2*r/(1+Q)+np.sqrt(a2/2*(1+1/Q)-(P*(1-e2)*z2)/(Q*(1+Q))-P*r2/2)
    U = np.sqrt((r-e2*r0)**2+z2)
    V = np.sqrt((r-e2*r0)**2+(1-e2)*z2)
    z0 = (b2*z)/(a*V)

    gdalt = U*(1-b2/(a*V))/1000.
    gdlat = np.arctan2(z+ep2*z0,r)*180./np.pi
    gdlon = 2*np.arctan2(r-x,y)*180./np.pi
    return gdlat, gdlon, gdalt

def geocentric_to_geodetic(gclat,gclon,gcalt):
    x, y, z = geocentric_to_cartesian(gclat,gclon,gcalt)
    gdlat, gdlon, gdalt = cartesian_to_geodetic(x,y,z)
    return gdlat, gdlon, gdalt

def geodetic_to_geocentric(gdlat,gdlon,gdalt):
    r, t, p = geodetic_to_spherical(gdlat,gdlon,gdalt)
    gclat, gclon, gcalt = spherical_to_geocentric(r,t,p)
    return gclat, gclon, gcalt

def spherical_to_geodetic(r,t,p):
    x, y, z = spherical_to_cartesian(r,t,p)
    gdlat, gdlon, gdalt = cartesian_to_geodetic(x,y,z)
    return gdlat, gdlon, gdalt

def geodetic_to_spherical(gdlat,gdlon,gdalt):
    lam_gd = gdlat*np.pi/180.
    h = gdalt*1000.
    phi = gdlon*np.pi/180.
    rho = Req/np.sqrt(1-e2*np.sin(lam_gd)**2)
    r = np.sqrt((rho+h)**2*np.cos(lam_gd)**2+(rho+h-e2*rho)**2*np.sin(lam_gd)**2)
    t = np.arccos((rho+h-e2*rho)*np.sin(lam_gd)/r)
    p = phi
    return r, t, p

def vector_spherical_to_cartesian(vr,vt,vp,r,t,p):
    vx = vr*np.sin(t)*np.cos(p)+vt*np.cos(t)*np.cos(p)-vp*np.sin(p)
    vy = vr*np.sin(t)*np.sin(p)+vt*np.cos(t)*np.sin(p)+vp*np.cos(p)
    vz = vr*np.cos(t)-vt*np.sin(t)
    return vx, vy, vz

def vector_cartesian_to_spherical(vx,vy,vz,x,y,z):
    r,t,p = cartesian_to_spherical(x,y,z)
    vr = vx*np.sin(t)*np.cos(p)+vy*np.sin(t)*np.sin(p)+vz*np.cos(t)
    vt = vx*np.cos(t)*np.cos(p)+vy*np.cos(t)*np.sin(p)-vz*np.sin(t)
    vp = -vx*np.sin(p)+vy*np.cos(p)
    return vr, vt, vp

def vector_spherical_to_geocentric(vr,vt,vp,*args):
    vnc = -vt
    vec = vp
    vuc = vr
    return vnc, vec, vuc

def vector_geocentric_to_spherical(vnc,vec,vuc,*args):
    vr = vuc
    vt = -vnc
    vp = vec
    return vr, vt, vp

def vector_geocentric_to_cartesian(vn,ve,vu,gclat,gclon,gcalt):
    r, t, p = geocentric_to_spherical(gclat,gclon,gcalt)
    vr, vt, vp = vector_geocentric_to_spherical(vn,ve,vu)
    vx, vy, vz = vector_spherical_to_cartesian(vr,vt,vp,r,t,p)
    return vx, vy, vz

def vector_cartesian_to_geocentric(vx,vy,vz,x,y,z):
    vr, vt, vp = vector_cartesian_to_spherical(vx,vy,vz,x,y,z)
    vn, ve, vu = vector_spherical_to_geocentric(vr,vt,vp)
    return vn, ve, vu


def vector_geocentric_to_geodetic(vnc,vec,vuc,gclat,gclon,gcalt):
    gdlat, gdlon, gdalt = geocentric_to_geodetic(gclat,gclon,gcalt)
    lam_gc = gclat*np.pi/180.
    lam_gd = gdlat*np.pi/180.
    b = lam_gd-lam_gc
    vud = vuc*np.cos(b)+vnc*np.sin(b)
    vnd = -vuc*np.sin(b)+vnc*np.cos(b)
    ved = vec
    return vnd, ved, vud

def vector_geodetic_to_geocentric(vnd,ved,vud,gdlat,gdlon,gdalt):
    gclat, gclon, gcalt = geodetic_to_geocentric(gdlat,gdlon,gdalt)
    lam_gc = gclat*np.pi/180.
    lam_gd = gdlat*np.pi/180.
    b = lam_gd-lam_gc
    vuc = vud*np.cos(-b)+vnd*np.sin(-b)
    vnc = -vud*np.sin(-b)+vnd*np.cos(-b)
    vec = ved
    return vnc, vec, vuc

def vector_spherical_to_geodetic(vr,vt,vp,r,t,p):
    gclat, gclon, gcalt = spherical_to_geocentric(r,t,p)
    vnc, vec, vuc = vector_spherical_to_geocentric(vr,vt,vp)
    vnd, ved, vud = vector_geocentric_to_geodetic(vnc,vec,vuc,gclat,gclon,gcalt)
    return vnd, ved, vud

def vector_geodetic_to_spherical(vnd,ved,vud,gdlat,gdlon,gdalt):
    vnc, vec, vuc = vector_geodetic_to_geocentric(vnd,ved,vud,gdlat,gdlon,gdalt)
    vr, vt, vp = vector_geocentric_to_spherical(vnc,vec,vuc)
    return vr, vt, vp

def vector_cartesian_to_geodetic(vx,vy,vz,x,y,z):
    gclat, gclon, gcalt = cartesian_to_geocentric(x,y,z)
    vnc, vec, vuc = vector_cartesian_to_geocentric(vx,vy,vz,x,y,z)
    vnd, ved, vud = vector_geocentric_to_geodetic(vnc,vec,vuc,gclat,gclon,gcalt)
    return vnd, ved, vud

def vector_geodetic_to_cartesian(vnd,ved,vud,gdlat,gdlon,gdalt):
    r, t, p = geodetic_to_spherical(gdlat,gdlon,gdalt)
    vr, vt, vp = vector_geodetic_to_spherical(vnd,ved,vud,gdlat,gdlon,gdalt)
    vx, vy, vz = vector_spherical_to_cartesian(vr,vt,vp,r,t,p)
    return vx, vy, vz


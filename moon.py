import numpy as np
import matplotlib.pyplot as plt
import math

def vec(a, b):
    return np.array([a,b],dtype=np.float64)
def vec3(a, b,c):
    return np.array([a,b,c],dtype=np.float64)
def normalize(a):
    return a/np.linalg.norm(a)

def rotationMatrix(axis, angle):
    cosA = np.cos(angle)
    sinA = np.sin(angle)
    print(angle.shape,'angleshape')
    matrixes = np.zeros([angle.shape[0],3,3])
    if axis == 0:

        return np.array([
            [1, 0, 0],
            [0, cosA, sinA],
            [0, -sinA, cosA]])
    if axis == 1:
        return np.array([
            [cosA, 0, -sinA],
            [0, 1, 0],
            [sinA, 0, cosA]])
    if axis == 2: 
        return np.array([
            [cosA, sinA, 0],
            [-sinA, cosA, 0],
            [0, 0, 1]])



class Orbit:
    def __init__(self,f1,f2,a,c,b,e,period,obliquity = 0,lon_adjustment = 0):
        self.f1 = f1
        self.f2 = f2
        self.a = a
        self.c = c
        self.b = b
        self.e = e
        self.period = period# in days
        self.obliquity = obliquity
        self.lon_adjustment = lon_adjustment
    
    def day_to_angle(self,day):
        return day/self.period*2*math.pi

    def goal_angle_to_orbital_pos(self, goal_angle):
        angle = goal_angle + 0
        M = goal_angle - self.e * np.sin(goal_angle - math.pi)
        goal_dif = M - goal_angle

        for n in range(10):
            angle += goal_dif
            M = angle - self.e * np.sin(angle - math.pi)
            goal_dif = goal_angle - M
        

        p = vec(np.cos(angle) * self.a, np.sin(angle) * self.b)
        return self.f1[:,np.newaxis] - p
    
    def planet_rotation_transform(self,day):
        day_angle = day * 366.25/365.25*2*np.pi
        day_matrix = rotationMatrix(2,day_angle)

        tilt_angle = np.deg2rad( self.obliquity )
        tilt_matrix = rotationMatrix(1, tilt_angle)

        angle_alignment_semi_major_to_tilt = -0.22363
        # to adjust for the fact that solstices do not align with periapsis and apoapsis

        day_tilt_to_elipse = rotationMatrix(2, angle_alignment_semi_major_to_tilt)

        full_matrix = np.dot(day_matrix, np.dot(tilt_matrix, day_tilt_to_elipse))
        return full_matrix
    
    def get_coords(self,point):
        p = normalize(point)
        lon = np.arctan2(p[1],p[0]) + self.lon_adjustment #grenich mean time adj
        lat = np.arctan2(p[2], np.sqrt(p[1] * p[1] + p[0] * p[0]))
        return vec(lon,lat)
    
    def get_transformed_point_at_days(self,points,days):
        angles = self.day_to_angle(days)
        points -= self.goal_angle_to_orbital_pos(angles)

        print(self.planet_rotation_transform(days).shape)

        



#sun orbit characteristics units of distance are R_earth
f1 = vec(-783.79 / 2, 0)
f2 = vec(783.79 / 2, 0)
a = 23455
c = abs(f1[0] - f2[0]) / 2
b = math.sqrt(a * a - c * c)
e = 0.0167086
obliquity = 23.44 #degree
period = 365.2425 #day
lon_adjustment=-0.2438 # radians this is just to start the day at grenich mean time
#may be a better way

EarthOrbit = Orbit(f1,f2,a,b,c,e,period,obliquity,lon_adjustment)

#angles = np.arange(100)*np.pi/100*2

#toplot = EarthOrbit.goal_angle_to_orbital_pos(angles)
points = np.zeros([2,100],dtype=np.float64)
days = np.arange(100)*365/100
toplot = EarthOrbit.get_transformed_point_at_days(points,days)

EarthOrbit.planet_rotation_transform(23.2)

# plt.plot(toplot[0,:],toplot[1,:])
# plt.plot(f1[0],f1[1],'ro')
# plt.plot(f2[0],f2[1],'ro')
# print(f1[0],f1[1],'ro')
# plt.show()

#moon orbit
perigree = 56.949
apogee = 63.561
a = (perigree+apogee)/2
c = (perigree-apogee)/2
f1 = vec(-c,0)
f2 = vec(c,0)
b = math.sqrt(a * a - c * c)
e = 0.0549006
obliquity = 6.68 #degree probably wont matter unless looking from moon or trying to get very accurate with how moon looks (wobble)
inclination = 5.14 #degree, matters a lot for where moon is in sky. angle to ecliptic plane
period = 27.321661 #days
#obliquity orbits at a rate of once every 18.6 years probably wont count that
# the foci/semimajor axis rotates about once every 8.85 years this changes the apparent diameter
# the inclination plane rotates by 360degrees once every 18.6 years with respect to the ecliptic

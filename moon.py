import numpy as np
import matplotlib.pyplot as plt
import math

def vec(a, b):
    return np.array([a,b],dtype=np.float64)
def vec3(a, b,c):
    return np.array([a,b,c],dtype=np.float64)
def normalize3(a):
    return a[:]/np.sqrt(a[:,0]*a[:,0]+a[:,1]*a[:,1]+a[:,2]*a[:,2])[:,np.newaxis]

def normalize(a):
    return a/np.linalg.norm(a)

def rotationMatrix(axis, angle):
    cosA = np.cos(angle)
    sinA = np.sin(angle)

    matrixes = np.zeros([angle.shape[0],3,3])
    if axis == 0:
        matrixes[:,1,1] = cosA
        matrixes[:,1,2] = sinA
        matrixes[:,2,1] = -sinA
        matrixes[:,2,2] = cosA
        matrixes[:,0,0] = 1
    if axis == 1:

        matrixes[:,0,0] = cosA
        matrixes[:,0,2] = -sinA
        matrixes[:,2,0] = sinA
        matrixes[:,2,2] = cosA
        matrixes[:,1,1] = 1
    if axis == 2: 

        matrixes[:,0,0] = cosA
        matrixes[:,0,1] = sinA
        matrixes[:,1,0] = -sinA
        matrixes[:,1,1] = cosA
        matrixes[:,2,2] = 1
    return matrixes
def rotationMatrix_single(axis, angle):
    cosA = np.cos(angle)
    sinA = np.sin(angle)
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
    def __init__(self,f1,f2,a,b,c,e,period,obliquity = 0,lon_adjustment = 0, inclination = 0,day_to_angle_adjustment = 0,inclination_to_ecliptic = 0):
        self.f1 = f1
        self.f2 = f2
        self.a = a
        self.c = c
        self.b = b
        self.e = e
        self.period = period# in days
        self.obliquity = obliquity
        self.lon_adjustment = lon_adjustment
        self.inclination = inclination#deg
        self.day_to_angle_adjustment = day_to_angle_adjustment
        self.inclination_to_ecliptic=inclination_to_ecliptic
    
    def day_to_angle(self,day,message = ''):
        # print('dayanglefunc',self.day_to_angle_adjustment)
        # print(message)
        return day/self.period*2*math.pi +self.day_to_angle_adjustment

    def angle_to_day(self,angle):
        # print('dayanglefunc',self.day_to_angle_adjustment)
        # print(message)
        return (angle)/(2*math.pi)*self.period

    def goal_angle_to_orbital_pos(self, goal_angle):
        angle = goal_angle + 0
        M = goal_angle - self.e * np.sin(goal_angle - math.pi)
        goal_dif = M - goal_angle

        for n in range(10):
            angle += goal_dif
            M = angle - self.e * np.sin(angle - math.pi)
            goal_dif = goal_angle - M
        

        p = np.stack([np.cos(angle) * self.a, np.sin(angle) * self.b,angle*0],axis=1)
        

        #     print(x)
        return self.f1 - p
    
    def goal_angle_to_orbital_pos_withinclination(self, goal_angle):
        p = self.goal_angle_to_orbital_pos( goal_angle)

        # inclination_matrix = rotationMatrix(1,np.deg2rad(self.inclination)+goal_angle*0)

    #     matrixes = np.zeros([goal_angle.shape[0],3,3])
    #     if axis == 0:
    #         matrixes[:,1,1] = cosA
    #         matrixes[:,1,2] = sinA
    #         matrixes[:,2,1] = -sinA
    #         matrixes[:,2,2] = cosA
    #         matrixes[:,0,0] = 1
    #     if axis == 1:

    #         matrixes[:,0,0] = cosA
    #         matrixes[:,0,2] = -sinA
    #         matrixes[:,2,0] = sinA
    #         matrixes[:,2,2] = cosA
    #         matrixes[:,1,1] = 1
    #         if axis == 1:
    #     ([
    #         [cosA, 0, -sinA],
    #         [0, 1, 0],
    #         [sinA, 0, cosA]])
    # if axis == 2: 
    #     return np.array([
    #         [cosA, sinA, 0],
    #         [-sinA, cosA, 0],
    #         [0, 0, 1]])

        #adjustment for Nodal precession
        day = self.angle_to_day(goal_angle)
        procession_angle = day / 6793 * 2 * np.pi

        inclination_to_ecliptic = self.inclination_to_ecliptic + procession_angle
        # inclination_to_ecliptic_matrix = rotationMatrix(2,inclination_to_ecliptic)
        
        # matrix = np.matmul(inclination_matrix,inclination_to_ecliptic_matrix)

        axis = np.stack([np.cos(inclination_to_ecliptic), 
                         np.sin(inclination_to_ecliptic),goal_angle*0],axis=1)
        #print(axis)
        a = np.deg2rad(self.inclination)
        #p = p*np.cos(a) + np.cross(axis,p,axis=1)*np.sin(a) + axis*np.dot(axis,p)*(1-np.cos(a))
        p1 =  p*np.cos(a)
        p2 = np.cross(axis,p,axis=1)*np.sin(a) 
        p3 = axis*(np.einsum('ij,ij->i',axis,p)*(1-np.cos(a)))[:,np.newaxis]

        p = p1+p2+p3
        #print(p)
        # print('matrix dims',matrix.shape)
        # p =  np.einsum('ijk,ik->ij', matrix, p)
        return p


    
    def planet_rotation_transform(self,day):
        day_angle = day * 366.25/365.25*2*np.pi
        day_matrix = rotationMatrix(2,day_angle)

        tilt_angle = np.deg2rad( self.obliquity )
        tilt_matrix = rotationMatrix_single(1, tilt_angle)

        angle_alignment_semi_major_to_tilt = -0.22363
        # to adjust for the fact that solstices do not align with periapsis and apoapsis

        day_tilt_to_elipse = rotationMatrix_single(2, angle_alignment_semi_major_to_tilt)

        tilt_and_el = np.dot(tilt_matrix, day_tilt_to_elipse)
        #day_matrix[1:] = 0
        full_matrix = np.dot(day_matrix, tilt_and_el)

        return full_matrix
    
    def get_coords(self,point):
        p = normalize3(point)
        lon = np.arctan2(p[:,1],p[:,0]) + self.lon_adjustment #grenich mean time adj
        lat = np.arctan2(p[:,2], np.sqrt(p[:,1] * p[:,1] + p[:,0] * p[:,0]))
        return vec(lon,lat)
    
    def get_coords_single(self,point):
        p = normalize(point)
        lon = np.arctan2(p[1],p[0]) + self.lon_adjustment #grenich mean time adj
        lat = np.arctan2(p[2], np.sqrt(p[1] * p[1] + p[0] * p[0]))
        return vec(lon,lat)
    
    def get_transformed_point_at_days(self,points,days):
        angles = self.day_to_angle(days)
        points -= self.goal_angle_to_orbital_pos(angles)
        

        matrixes = self.planet_rotation_transform(days)

        #p = np.dot(matrix,points)
        #p = np.tensordot(matrix, points, axes=([1], [1]))

        #p = points[:, 0:1] * matrixes[:, :, 0] + points[:, 1:2] * matrixes[:, :, 1] + points[:, 2:3] * matrixes[:, :, 2]
        p =  np.einsum('ijk,ik->ij', matrixes, points)
        # p1 = np.sum(points[:] * matrixes[:, 0], axis = 1)
        # p2 = np.sum(points[:] * matrixes[:,1], axis = 1)
        # p3 = np.sum(points[:] * matrixes[:,2], axis = 1)
        # p =np.stack( [p1,p2,p3], axis = 1) 
        #p = normalize3(p)


        return normalize3(p)


        



#sun orbit characteristics units of distance are R_earth
f1 = vec3(-783.79 / 2, 0,0)
f2 = vec3(783.79 / 2, 0,0)
a = 23455
c = abs(f1[0] - f2[0]) / 2
b = math.sqrt(a * a - c * c)
e = 0.0167086
obliquity = 23.44 #degree
period = 365.2425 #day
lon_adjustment=-0.2438 # radians this is just to start the day at grenich mean time
#may be a better way
dayadjustment = 182/365.25*np.pi*2

EarthOrbit = Orbit(f1,f2,a,b,c,e,period,obliquity,lon_adjustment,day_to_angle_adjustment=dayadjustment)

#angles = np.arange(100)*np.pi/100*2


n = 365*12
points = np.zeros([1,3],dtype=np.float64)
days = np.arange(n)/(365.25-12.37)*365.25
dayse = np.array([255 + 0.70898])
earthangles = EarthOrbit.day_to_angle(dayse)
earthloc = EarthOrbit.goal_angle_to_orbital_pos(earthangles)
toplot = EarthOrbit.get_transformed_point_at_days(points,dayse)

print(EarthOrbit.get_coords(toplot),'coords')

# #EarthOrbit.planet_rotation_transform(23.2)
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(toplot[:,0],toplot[:,1],toplot[:,2])
#plt.scatter(toplot[0,:],toplot[1,:],toplot[3,:])
# plt.plot(f1[0],f1[1],'ro')
# plt.plot(f2[0],f2[1],'ro')

#plt.show()

#moon orbit
perigree = 56.949
apogee = 63.561
a = (perigree+apogee)/2
c = (perigree-apogee)/2
f1 = vec3(-c,0,0)
f2 = vec3(c,0,0)
b = math.sqrt(a * a - c * c)
e = 0.0549006
obliquity = 6.68 #degree probably wont matter unless looking from moon or trying to get very accurate with how moon looks (wobble)
inclination = 5.14 #degree, matters a lot for where moon is in sky. angle to ecliptic plane
period = 27.321661 #days
#obliquity orbits at a rate of once every 18.6 years probably wont count that
# the foci/semimajor axis rotates about once every 8.85 years this changes the apparent diameter
# the inclination plane rotates by 360degrees once every 18.6 years with respect to the ecliptic

day_to_angle_adjustment = 2.638937829015426


inclination_to_ecliptic = 4.335397861953914

MoonOrbit = Orbit(f1,f2,a,b,c,e,period,obliquity,inclination=inclination,
                  day_to_angle_adjustment=day_to_angle_adjustment,
                  inclination_to_ecliptic=inclination_to_ecliptic)

# moonangle = MoonOrbit.day_to_angle(days)
# moonloc = MoonOrbit.goal_angle_to_orbital_pos_withinclination(moonangle)+earthloc

# toplot2 = EarthOrbit.get_transformed_point_at_days(moonloc,days)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(toplot[:,0],toplot[:,1],toplot[:,2])
ax.axes.set_xlim3d(left=-1, right=1) 
ax.axes.set_ylim3d(bottom=-1, top=1) 
ax.axes.set_zlim3d(bottom=-1, top=1 ) 
# num_items = toplot2.shape[0]
# colors = np.arange(num_items) // 29
# ax.scatter(toplot2[:,0],toplot2[:,1],toplot2[:,2],c=colors,cmap='hsv')

def get_moon_loc(days):
    #points = np.zeros([days.shape[0],3],dtype=np.float64)
    earthangles = EarthOrbit.day_to_angle(days)
    earthloc = EarthOrbit.goal_angle_to_orbital_pos(earthangles)
    moonangle = MoonOrbit.day_to_angle(days,'moon')
    moonloc = MoonOrbit.goal_angle_to_orbital_pos_withinclination(moonangle)+earthloc
    toplot2 = EarthOrbit.get_transformed_point_at_days(moonloc,days)
    return toplot2

def get_sun_loc(days):
    #points = np.zeros([days.shape[0],3],dtype=np.float64)
    earthangles = EarthOrbit.day_to_angle(days)
    earthloc = EarthOrbit.goal_angle_to_orbital_pos(earthangles)
    
    return earthloc

def plotmoon3(days):
    toplot2 = get_moon_loc(days)
    ax.scatter(toplot2[:,0],toplot2[:,1],toplot2[:,2],c=days,cmap='hsv')
    print("moon coords",EarthOrbit.get_coords(toplot2))


#d1 = np.arange(13)*27.321661
#plotmoon3(days)
today = np.array([254.91667])
#plotmoon3(today)

plotmoon3(days)

plt.show()
print(get_moon_loc(np.array([97])),get_sun_loc(np.array([97])))

##sl = get_sun_loc(dayse)
#print(MoonOrbit.get_coords(sl))

best = [1000000,10000000]
closest = 10000000000000000

goal = vec(0.31881,-0.49626)

for x in range(0,200):
    for y in range(-100,100):
        MoonOrbit=0
        MoonOrbit = Orbit(f1,f2,a,b,c,e,period,obliquity,inclination=inclination,
        day_to_angle_adjustment = x/100*np.pi*2,inclination_to_ecliptic = y/100*np.pi*2)
        #print(MoonOrbit.inclination_to_ecliptic)
        A = x/100*np.pi*2
        B = y/100*np.pi*2
        sl = get_moon_loc(dayse)
        #print(sl)
        c = MoonOrbit.get_coords(sl)
        dist = np.linalg.norm(c-goal)
        if dist <= closest:
            print("best")
            closest = dist + 0
            print(c,A,B,dist)

        #print(c)

        
        
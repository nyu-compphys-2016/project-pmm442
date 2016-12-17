import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll

d = 3.57 #mu m
#d = 12.5 #mu m gives circle 
a = 1.0 #mu m
L = 10.0 #mu m
v0 = 40.0 #mu m/s
#I0 = 0.1 #W
kb = 1.38e-11 #kg * mu m^2 * s^-2 * K^-1
T = 294 #K room temperature
visw = 1.002e-9 # viscosity of water kg * m * s^-1
Diff = (1/(6*np.pi*visw*a))*kb*T
#gom = 100.0
c = 3e8 # m/s
#alphar = 3 #mu m^3/s/mW
#mutheta = 

def intensity(r,I0=1.0):
	sigma = 7.3 #mu m
	return I0*np.exp(-(r*r)/(2*sigma*sigma))

def evolve_sigma(h,noise):
 	xw, yw, phiw = np.random.normal(size = 3, loc=0.0, scale=np.sqrt(Diff*h)) *noise
 	phiw *= np.sqrt(3/(4*a*a))
 	return np.array([xw, yw, phiw])

def evolve_mu(vec):
	x = vec[0]
	y = vec[1]
	phi = vec[2]
	r = np.sqrt(x*x + y*y)
	v = v0*intensity(r)
	xdot = v*np.cos(phi)
	ydot = v*np.sin(phi)
	phidot = (x*ydot - y*xdot)/(d*np.sqrt(r*r + L*L))
	return np.array([xdot,ydot,phidot])

def euler_mar(vec,tf,h,noise):
	xpoints = []
	ypoints = []
	phipoints = []
	t = 0.0
	while t < tf:
		if t + h >tf:
			h = tf - t
		xpoints.append(vec[0])
		ypoints.append(vec[1])
		phipoints.append(vec[2])
		dvec = evolve_mu(vec)
		vec = vec + h*dvec + evolve_sigma(h,noise)
		t = t + h
	xpoints.append(vec[0])
	ypoints.append(vec[1])
	phipoints.append(vec[2])
	return xpoints, ypoints, phipoints, t 
	
def withoutm(phi0, tf):
	x0 = 15.0
	y0 = 0.0
	#phi0 = np.pi/2
	
	N = 1000000.0
	I0 = 1.0
	#tf = 60.0				#s
	h = (tf)/N
	
	#vec = np.array([x0,y0,phi0])
	#xnonoise,ynonoise,phinonoise,tnonoise=euler_mar(vec,tf,h,0.0)
	
	vec = np.array([x0,y0,phi0])
	xnoise,ynoise,phinoise,tnoise=euler_mar(vec,tf,h,1.0)
	
	#lc1 = colorline(xnonoise,ynonoise,cmap='inferno')
	#plt.xlabel(r'$x$ [$\mu$m]')
	#plt.ylabel(r'$y$ [$\mu$m]')
	#cb1 = plt.colorbar(lc1)
	#cb1.set_label(r'$t$ [s]')
	#plt.grid(True)
	#plt.title(r'No thermal noise')
	#plt.axis('equal')
	
	#plt.show()
	
	lc2 = colorline(xnoise,ynoise,cmap='inferno')
	plt.axis('equal')
	plt.xlabel(r'$x$ [$\mu$m]')
	plt.ylabel(r'$y$ [$\mu$m]')
	cb2 = plt.colorbar(lc2)
	cb2.set_label(r'$t$ [s]')
	plt.grid(True)
	plt.title(r'Stochastic thermal noise')
	plt.axis('equal')
	
	plt.show()

def mevolve_mu(vec,gom):
	x = vec[0]
	y = vec[1]
	phi = vec[2]
	vx = vec[3]
	vy = vec[4]
	r = np.sqrt(x*x + y*y)
	v = v0*intensity(r)
	xdot = v*np.cos(phi)
	ydot = v*np.sin(phi)
	vxdot = gom*(xdot - vx)
	vydot = gom*(ydot - vy)
	phidot = (x*ydot - y*xdot)/(d*np.sqrt(r*r + L*L))
	return np.array([vx,vy,phidot,vxdot,vydot])

def meuler_mar(vec,tf,h,gom,noise):
	xpoints = []
	ypoints = []
	phipoints = []
	t = 0.0
	while t < tf:
		if t + h >tf:
			h = tf - t
		xpoints.append(vec[0])
		ypoints.append(vec[1])
		phipoints.append(vec[2])
		dvec = mevolve_mu(vec,gom)
		vec[0:2] = vec[0:2] + h*dvec[0:2]
		vec[2:5] = vec[2:5] + h*dvec[2:5] + evolve_sigma(h,noise)
		t = t + h
	xpoints.append(vec[0])
	ypoints.append(vec[1])
	phipoints.append(vec[2])
	return xpoints, ypoints, phipoints, t 
	
def withm(phi0, tf):
	x0 = 15.0
	y0 = 0.0
	vx = 0.0
	vy = 40.0
	#phi0 = np.pi/2
	#noise = 1.0
	
	N = 1000000.0
	#tf = 60.0				#s
	h = (tf)/N
	
	for j in range(6):
		if j == 0:
			noise = 0.0
		else:
			noise = 1.0
		for i in range(4):
			gom = 10**i
		
			vec = np.array([x0,y0,phi0,vx,vy])
			xpoints,ypoints,phipoints,t=meuler_mar(vec,tf,h,gom,noise)
	
			lc = colorline(xpoints,ypoints,cmap='inferno',tf=tf)
			plt.axis('equal')
			plt.xlabel(r'$x$ [$\mu$m]')
			plt.ylabel(r'$y$ [$\mu$m]')
			cb = plt.colorbar(lc)
			cb.set_label(r'$t$ [s]')
			plt.grid(True)
			plt.axis('equal')
			if noise == 0.0:
				plt.title(r'No noise with $\frac{\gamma}{m}$=%.1f'%gom)
				plt.savefig('/Users/paulmcnulty/Desktop/PHYSICS/comp_phys/massnonoise%03dlong.png'%i)
			else:
				plt.title(r'Stochastic thermal noise with $\frac{\gamma}{m}$=%.1f'%gom)
				plt.savefig('/Users/paulmcnulty/Desktop/PHYSICS/comp_phys/massnoise%03dlong%03d.png'%(i,j))
			plt.close()
	#withi(phi0,tf)

def ievolve_mu(vec,I0):
	x = vec[0]
	y = vec[1]
	phi = vec[2]
	r = np.sqrt(x*x + y*y)
	I = intensity(r,I0)
	v = v0*I/I0
	xdot = v*np.cos(phi) + np.pi*I*a*a*np.cos(np.arctan2(y,x))
	ydot = v*np.sin(phi) + np.pi*I*a*a*np.sin(np.arctan2(y,x))
	phidot = (x*ydot - y*xdot)/(d*np.sqrt(r*r + L*L))
	return np.array([xdot,ydot,phidot])

def ieuler_mar(vec,tf,h,I0,noise):
	xpoints = []
	ypoints = []
	phipoints = []
	t = 0.0
	while t < tf:
		if t + h >tf:
			h = tf - t
		xpoints.append(vec[0])
		ypoints.append(vec[1])
		phipoints.append(vec[2])
		dvec = ievolve_mu(vec,I0)
		vec = vec + h*dvec + evolve_sigma(h,noise)
		t = t + h
	xpoints.append(vec[0])
	ypoints.append(vec[1])
	phipoints.append(vec[2])
	return xpoints, ypoints, phipoints, t 
	
def withi(phi0, tf):
	x0 = 15.0
	y0 = 0.0
	#phi0 = np.pi/2
	noise = 1.0
	
	N = 1000000.0
	#I0 = 1.0e8
	#tf = 60.0				#s
	h = (tf)/N
	
	I0s = np.arange(0.1,2.0,0.5)
	for j in range(6):
		if j == 0:
			noise = 0.0
		else:
			noise = 1.0
		for I0 in I0s:
			vec = np.array([x0,y0,phi0])
			xpoints,ypoints,phipoints,t=ieuler_mar(vec,tf,h,I0,noise)
	
			lc = colorline(xpoints,ypoints,cmap='inferno',tf=tf)
			plt.axis('equal')
			plt.xlabel(r'$x$ [$\mu$m]')
			plt.ylabel(r'$y$ [$\mu$m]')
			cb = plt.colorbar(lc)
			cb.set_label(r'$t$ [s]')
			plt.grid(True)
			plt.axis('equal')
			if noise == 0.0:
				plt.title(r'No noise with scaled $I_0$=%.2f'%I0)
				plt.savefig('/Users/paulmcnulty/Desktop/PHYSICS/comp_phys/intennonoise%.3flong.png'%I0)
			else:
				plt.title(r'Stochastic thermal noise with scaled $I_0$=%.2f'%I0)
				plt.savefig('/Users/paulmcnulty/Desktop/PHYSICS/comp_phys/intennoise%.3flong%d.png'%(I0,j))
			plt.close()
		
def colorline(x,y,z=None,cmap='copper',tf=60.0,linewidth=2,alpha=1.0):
	norm=plt.Normalize(0.0,tf)
	
	if z is None:
		z = np.linspace(0.0,tf,len(x))
		
	z = np.asarray(z)
	
	segments = make_segments(x,y)
	lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
	
	ax = plt.gca()
	ax.add_collection(lc)
	
	return lc
	
def make_segments(x,y):
	points = np.array([x,y]).T.reshape(-1,1,2)
	segments = np.concatenate([points[:-1],points[1:]],axis=1)
	return segments
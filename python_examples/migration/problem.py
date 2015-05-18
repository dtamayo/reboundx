import rebound
import reboundxf
import math

#gamma=1.5
#Rc=80.
#podot=0.01
rebound.G = 4*math.pi**2
starmass = 0.55
mass = 5.5e-4

rebound.add(m=starmass)

#alpha = G*starmass/Rc/Rc*podot/(1.-gamma/2.)
a0s = [46.,77.3]

for a in a0s:
    #vcirc = math.sqrt(G*starmass/a + alpha*Rc*(Rc/a)**(gamma-1.))*1.05
    rebound.add(m=mass,a=a,e=0.)

rebound.additional_forces = reboundxf.forces()
reboundxf.set_e_damping([1.e2,1.e2,1.e2])
#reboundxf.set_peri_precession(gamma,Rc,podot)
reboundxf.set_migration([0.,-1.e6,0.])

rebound.integrator = "ias15"
rebound.integrate(1.)

'''
particles = rebound.particles_get()

tmax = 1000000.

last_t = -1e6
outputdelta = 0.1

a = []
e = []
i = []
po = []
Omega = []
omega = []
t = []
x = []
x2 = []
r = []
r2 = []

while rebound.get_t() < tmax:
    _t = rebound.get_t()
    if _t - last_t > outputdelta:
        #o = pytools.p2orbit(particles[1], particles[0])
        t.append(_t)
        x.append(particles[1].x)
        x2.append(particles[2].x)
        r.append(math.sqrt(particles[1].x**2 + particles[1].y**2 + particles[1].z**2))
        r2.append(math.sqrt(particles[2].x**2 + particles[2].y**2 + particles[2].z**2))
        a.append(o.a)
        e.append(o.e)
        i.append(o.inc)
        po.append(o.Omega + o.omega)
        last_t = _t

    rebound.step()
    
print(x[10000:10100])
print(np.mean(x[:100]))
print(np.mean(x[:1000]))
print(np.mean(x[:10000]))
wind = 100000
rmean = np.convolve(r, np.ones((wind,))/wind, mode='valid')
xmean = np.convolve(x, np.ones((wind,))/wind, mode='valid')
tmean = np.convolve(t, np.ones((wind,))/wind, mode='valid')


tsmoothed = [np.mean(i) for i in np.array_split(np.asarray(t),1000.)]
rmax = [max(i) for i in np.array_split(np.asarray(r),1000.)]
rmin = [min(i) for i in np.array_split(np.asarray(r),1000.)]

a = [(rmax[i] + rmin[i])/2 for i in range(len(rmax))]
e = [(rmax[i] - rmin[i])/2/(rmax[i] + rmin[i]) for i in range(len(rmax))]

rmax2 = [max(i) for i in np.array_split(np.asarray(r2),1000.)]
rmin2 = [min(i) for i in np.array_split(np.asarray(r2),1000.)]

a2 = [(rmax2[i] + rmin2[i])/2 for i in range(len(rmax2))]
e2 = [(rmax2[i] - rmin2[i])/2/(rmax2[i] + rmin2[i]) for i in range(len(rmax2))]

fig, ax = plt.subplots(4)

#ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[2].set_yscale('log')

ax[0].plot(t[::100],a[::100])
ax[1].plot(t[::100],e[::100])
#ax[2].plot(t,i)
ax[3].plot(t,po, '.')


fig,ax = plt.subplots(2)
ax[1].set_yscale('log')
ax[0].plot(tsmoothed,a)
ax[0].plot(tsmoothed,a2)
ax[1].plot(tsmoothed,e)
ax[1].plot(tsmoothed,e2)

#pytools.plot_freq_spectrum(np.asarray(t[len(t)//2:]), np.asarray(x[len(t)//2:]),2.e-4, 1.e-3,nfreq=2000,log=False)
plt.show()
  '''      

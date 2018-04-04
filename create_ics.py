import numpy as np

from amuse.lab import *

from parameters import Parameters
### First some basics

def semimajoraxis(
        Mother,
        Mself,
        period,
        ):
    a3 = (
            (period / (2*np.pi))**2 * \
            constants.G * (Mother + Mself) \
            )
    return a3**(1/3.)# * (1 - Mself/Mother)

#def kepler_orbit(
#        orbiters, 
#        primary, 
#        a_over_r        = 1.,
#        direction       = "prograde",
#        eccentricity    = 0,
#        ):
#    """
#    Calculates velocity of the orbiters to be in Kepler orbit around
#    the primary, with the set eccentricity.
#    """
#
#    if eccentricity != 0:
#        print "Nonzero eccentricity not yet supported"
#        return -1
#    
#    orbiters.r  = orbiters.position.lengths()
#
#    a           = a_over_r * orbiters.r
#    mu          = constants.G * primary.mass
#
#    orbiters.vx +=  (np.sin(orbiters.phi) * (mu * (2./orbiters.r - 1./a))**0.5)
#    orbiters.vy += -(np.cos(orbiters.phi) * (mu * (2./orbiters.r - 1./a))**0.5)
#    orbiters.vz += 0.0 | units.kms
#
#    orbiters.initial_orbital_period = 2 * np.pi * (a**3 / mu)**0.5
#
#    if direction == "prograde":
#        #logger.info("Prograde orbit of XXXX")
#        print "Prograde orbit of XXXX"
#    elif direction == "retrograde":
#        print "Retrograde orbit of XXXX"
#        #logger.info("Retrograde orbit of XXXX")
#        orbiters.velocity *= -1
#
#    return 1


def xyz_from_cylindrical(
        orbiters,
        ):
    orbiters.x = (orbiters.r * np.cos(orbiters.phi))
    orbiters.y = (orbiters.r * np.sin(orbiters.phi))
    return

def kepler_orbit(
        orbiters, 
        centre, 
        direction   = "prograde",
        ):
    orbiters.r = (orbiters.position - centre.position).lengths()
    a = orbiters.r / (1+orbiters.eccentricity)
    mu = constants.G * (centre.mass+orbiters.mass)
    if direction == "retrograde":
        d = -1
    else:
        d = 1
    orbiters.vx +=  d * (np.sin(orbiters.phi)*(mu * (2./orbiters.r - 1./a))**0.5)
    orbiters.vy += -d * (np.cos(orbiters.phi)*(mu * (2./orbiters.r - 1./a))**0.5)
    orbiters.vz +=  d * 0|units.kms
    orbiters.initial_orbital_period = 2 * np.pi * (a**3 / mu)**0.5
    return

def kepler_periastron_velocity(
        star,
        planet,
        epsilon = 0.0,
        a       = 1.0 | units.AU,
        ): 
    mu      = constants.G * (star.mass+planet.mass)
    v_per   = np.sqrt( ((1+epsilon)*mu)/((1-epsilon)*a))
    return v_per

def kepler_apastron_velocity(
        star,
        planet,
        epsilon = 0.0,
        a       = 1.0 | units.AU,
        ): 
    mu      = constants.G * (star.mass+planet.mass)
    v_apo   = np.sqrt( ((1-epsilon)*mu)/((1+epsilon)*a)) 
    return v_apo


### Now, some models

def new_planetary_system(
        p,
        ):
    if p.model_type == "solarsystem":
        from amuse.ext.solarsystem import new_solar_system
        return new_solar_system()[0:9]

    particles = Particles()

    # First, we create the massive particles in the system. 
    # We also set their parameters, as given.

    for i in range(p.stars):
        star                = Particle()
        star.type           = "star"
        star.number         = i
        star.mass           = p.stars_mass[i]
        star.radius         = p.stars_radius[i]
        star.eccentricity   = p.stars_eccentricity[i]
        star.semimajoraxis  = p.stars_semimajoraxis[i]
        particles.add_particle(star)

    for i in range(p.planets):
        planet                  = Particle()
        planet.type             = "planet"
        planet.number           = i
        planet.mass             = p.planets_mass[i]
        planet.radius           = p.planets_radius[i]
        planet.hoststar         = p.planets_host[i]
        planet.eccentricity     = p.planets_eccentricity[i]
        planet.semimajoraxis    = p.planets_semimajoraxis[i]
        particles.add_particle(planet)

    for i in range(p.moons):
        moon                = Particle()
        moon.type           = "moon"
        moon.number         = i
        moon.mass           = p.moons_mass[i]
        moon.radius         = p.moons_radius[i]
        moon.hostplanet     = p.moons_host[i]
        moon.eccentricity   = p.moons_eccentricity[i]
        moon.semimajoraxis  = p.moons_semimajoraxis[i]
        particles.add_particle(moon)

    # Now, we get the subsets
    stars      = particles.select(lambda t: t == "star", ["type"])
    planets    = particles.select(lambda t: t == "planet", ["type"])
    moons      = particles.select(lambda t: t == "moon", ["type"])

    # Next, we position them correctly and give their velocities
    for star in stars:
        star.position   = p.stars_position[i]
        star.velocity   = p.stars_velocity[i]

    for planet in planets:
        star    = stars.select(lambda n: n == planet.hoststar, ["number"])[0]
        planet.eccentricity = p.planets_eccentricity[i]
        planet.position     = np.array(
                [
                    (1+planet.eccentricity) * p.planets_semimajoraxis[i].value_in(units.AU), 
                    0, 
                    0,
                    ],
                ) | units.AU
        planet.position    += star.position
        planet.phi          = 0 ##fixme
        kepler_orbit(
                planet,
                star,
                )
        planet.velocity    += star.velocity

    for moon in moons:
        planet  = planets.select(lambda n: n == moon.hostplanet, ["number"])[0]
        moon.position   = np.array(
                [
                    p.moons_semimajoraxis[i].value_in(units.AU), 
                    0, 
                    0,
                    ],
                ) | units.AU
        moon.position  += planet.position
        moon.velocity  += planet.velocity
        kepler_orbit(
                moon,
                planet,
                direction   = p.moons_direction[i],
                #eccentricity    = p.moons_eccentricity[i],
                )

    # Finally, we add the discs
    for i in range(p.discs):
        hosttype   = p.discs_host[i][0]
        hostnumber = p.discs_host[i][1]
        if hosttype == "star":
            host    = stars.select(lambda n: n == hostnumber, ["number"])
        elif hosttype == "planet":
            host    = planets.select(lambda n: n == hostnumber, ["number"])
        else:
            logger.error("No such host: %s %i"%(hosttype, hostnumber))
                    
        disc                = new_disc_model(
                i,
                host,
                p,
                )
        disc.type           = "disc"
        disc.number         = i
        disc.hosttype       = p.discs_host[i][0]
        disc.hostnumber     = p.discs_host[i][1]
        particles.add_particles(disc)

    disc       = particles.select(lambda t: t == "disc", ["type"])

    return particles

def new_disc_model(
        i,
        host,
        p,
        ):
    if p.discs_model[i]     == "random":
        particles   = new_random_disc_model(
                i,
                p,
                )
    elif p.discs_model[i]   == "powerlaw":
        particles   = new_powerlaw_disc_model(
                i,
                p,
                #nrings      = nrings,
                #startring   = int(nrings * ring_min_dist/ring_max_dist),
                #scale       = ring_max_dist,
                )
    else:
        logger.error("No such disc model type: %s"%p.model_type)

    particles.position += host.position
    particles.velocity += host.velocity
    kepler_orbit(
        particles,
        host, 
        direction   = p.discs_direction[i],
        )
    rnd_vx_factor = 1 + ((np.random.random(len(particles))-0.5) * 0.05)
    rnd_vy_factor = 1 + ((np.random.random(len(particles))-0.5) * 0.05)

    #particles.vx *= rnd_vx_factor
    #particles.vy *= rnd_vy_factor

    return particles

def new_powerlaw_disc_model(
        i,
        p,
        alpha       = 0,
        mode        = "scattered",
        ):
        #):
    """
    Model for a largely uniform (or powerlaw) particle distribution in
    a disk, based on the work of Satoko Yamamoto
    """
    
    nrings      = p.discs_nrings[i]
    startring   = int(nrings*(p.discs_rmin[i]/p.discs_rmax[i]))
    scale       = p.discs_rmax[i]
    varp        = int(np.floor(np.pi * (alpha + 2)))

    particles = Particles()

    for ring in range(startring, nrings):
        n_in_ring   = 1 + ring * varp
        angle       = 2. * np.pi / n_in_ring
        
        this_ring = (Particles(n_in_ring))
        
        this_ring.phi    = angle * np.array(list(range(n_in_ring)))
        this_ring.ring   = ring
        particles.add_particles(this_ring)
   
    particles.eccentricity  = 0
    particles.r         = scale * np.array( particles.ring / float(nrings) )
    if mode == "scattered":
        particles.r        += (np.random.random(len(particles)) - 0.5) * (1./nrings) * scale 
    particles.mass      = 0|units.kg
    particles.radius    = 10|units.km
    particles.rstart    = particles.r
    particles.z         = (np.random.random(len(particles)) - 0.5) * p.discs_zrange[i]
    particles.velocity  = np.array([0.,0.,0.]) | units.kms
    xyz_from_cylindrical(particles)
    #particles.x = particles.r * -1 * np.cos(particles.phi)
    #particles.y = particles.r * np.sin(particles.phi)
    #particles.z = 0|units.AU


    return particles

def _new_planetary_system(
        solar_system            = False,
        star_mass               = 1.0 | units.MSun,
        star_radius             = 1.0 | units.RSun,
        number_of_planets       = 1,
        orbital_eccentricities  = 0,
        planet_masses           = 1.0 | units.MJupiter,
        planet_radii            = 1.0 | units.RJupiter,
        planet_semi_major_axes  = 4.0 | units.AU,
        ):
    if solar_system:
        return new_solar_system()[0:9]
    particles       = Particles(number_of_planets + 1)
    star            = particles[0]
    star.position   = quantity([0,0,0])|units.AU
    star.velocity   = quantity([0,0,0])|units.kms

    planets         = particles[1:]
    for i in range(len(planets)):
        planet          = planets[i]
        planet.mass     = planet_masses[i]
        planet.radius   = planet_radii[i]
        planet.position = quantity[planet_semi_major_axes[i],0,0]
        kepler_orbit(
                planet, 
                star,
                #eccentricity    = orbital_eccentricities[i],
                )

    return particles

def new_random_disc_model(
        i,
        p,
        ):

    discsize = p.discs_rmax[i] - p.discs_rmin[i]

    particles           = Particles(p.discs_n[i])
    particles.phi       = np.random.random(len(particles)) * 2 * np.pi
    particles.r         = np.random.random(len(particles)) * discsize + p.discs_rmin[i]
    particles.density   = p.discs_density[i]#3 | units.g * units.cm**-3
    particles.mass      = p.discs_mass[i] / len(particles)
    particles.radius    = 0.1 | units.km 
    particles.volume    = (particles.mass / particles.density)
    particles.rstart    = particles.r
    particles.eccentricity  = 0
    particles.z         = (np.random.random(len(particles)) - 0.5) * p.discs_zrange[i]
    particles.velocity  = np.array([0.,0.,0.]) | units.kms
    xyz_from_cylindrical(particles)
    #particles.radius    = 100|units.m#1|units.RJupiter#00| units.km#RJupiter

    ## change orientation of the ring disk with respect to the planetary orbit
    #x   = np.sin(theta) * particles.x
    #z   = np.cos(theta) * particles.x
    #vx  = np.sin(theta) * particles.vx
    #vz  = np.cos(theta) * particles.vx
    #particles.x     = x
    #particles.z     = z
    #particles.vx    = vx
    #particles.vz    = vz

    return particles

def new_random_uniform_disc_model(
        i,
        p,
        ):

    discsize = p.discs_rmax[i] - p.discs_rmin[i]

    particles           = Particles(p.discs_n[i])
    particles.phi       = np.random.random(len(particles)) * 2 * np.pi
    particles.r         = np.random.random(len(particles)) * discsize + p.discs_rmin[i]
    particles.density   = p.discs_density[i]#3 | units.g * units.cm**-3
    particles.mass      = p.discs_mass[i] / len(particles)
    particles.radius    = 0.1 | units.km 
    particles.volume    = (particles.mass / particles.density)
    particles.rstart    = particles.r
    particles.eccentricity  = 0
    particles.z         = (np.random.random(len(particles)) - 0.5) * p.discs_zrange[i]
    particles.velocity  = np.array([0.,0.,0.]) | units.kms
    xyz_from_cylindrical(particles)
    #particles.radius    = 100|units.m#1|units.RJupiter#00| units.km#RJupiter

    ## change orientation of the ring disk with respect to the planetary orbit
    #x   = np.sin(theta) * particles.x
    #z   = np.cos(theta) * particles.x
    #vx  = np.sin(theta) * particles.vx
    #vz  = np.cos(theta) * particles.vx
    #particles.x     = x
    #particles.z     = z
    #particles.vx    = vx
    #particles.vz    = vz

    return particles

if __name__ in "__main__":


    #nrings          = 160
    #relative_width  = 0.1
    #startring       = int((1./relative_width) * nrings)

    for model in ["A","B","C"]:
        for mass in [20,40,60,80,90,100]:
            for direction in ["P","R"]:
                p   = Parameters(model_type="%s%s%i"%(model,direction,mass))
                print(p.model_type, p.model_name) 

                particles   = new_planetary_system(
                        p,
                        #nrings      = nrings,
                        #startring   = int(nrings * 0.2),
                        #scale       = (0.06) | units.AU,
                        )
                print(len(particles), p.discs_rmax[0], p.planets_semimajoraxis[0])
                particles.collection_attributes.model_name = p.model_name
                particles.collection_attributes.model_type = p.model_type
                write_set_to_file(particles, "model-%s.hdf5"%(p.model_name), "amuse")

import os,sys
import numpy as np
import time as realtime

from amuse.lab import *
import logging

#from planet.tools import cylindrical_from_xyz

analysislogger = logging.getLogger("analysis")

set_printing_strategy(
        "custom", 
        preferred_units = [
            units.MJupiter, 
            units.AU, 
            units.yr, 
            units.kms,
            ], 
        precision   = 5,
        )
def cleanup(
        grav, 
        particles, 
        converter,
        system_max_distance = 1.|units.AU,
        box_max_distance    = 1000|nbody_system.length,
        ):
    return particles[1:]

def run_analysis(
        system,
        p,
        RH                  = 1|units.AU,
        boundradii          = False,
        ):
    
    analysislogger.debug(
            "#Running analyse at time %s"%(
                system.model_time,
                ),
            )

    planet_distance_to_star = (system.planets[0].position - system.stars[0].position).length()
    planet_orbital_phase = np.arctan2(
            (system.planets[0].y - system.stars[0].y).value_in(units.AU),
            (system.planets[0].x - system.stars[0].x).value_in(units.AU),
            )

    ### number of bound particles?!
    # calculate Kepler properties for all ring particles
    # if e > 1: unbound, else bound

    #r   = analyse(
    #        system.disc, 
    #        system.planets[0], 
    #        len(system.ring), 
    #        RH          = RH, 
    #        boundradii  = boundradii,
    #        )

    

    analysislogger.info(
            "%s %s %s %s %s %s %s %s %s %s"%(
                system.model_time,
                system.gravity.model_time,
                planet_orbital_phase,
                len(system.particles),
                (system.planets[0].x-system.stars[0].x),
                (system.planets[0].y-system.stars[0].y),
                (system.planets[0].vx-system.stars[0].vx),
                (system.planets[0].vy-system.stars[0].vy),
                system.disc.y.max()-system.disc.y.min(),
                (system.disc.y.max()-system.disc.y.min())/(system.planets[0].vy-system.stars[0].vy),
                )
            )

    #if r50 > ring_max_distance: ##FIXME
    #    exit_dissolving()

    return 

def __analyse(
        ring,
        planet,
        ring_shells,
        ):
    ring.r = (ring.position - planet.position).lengths()

    rmin = []
    rmax = []
    ravg = []
    rstd = []
    rmed = []

    for shell in ring_shells:
        shell = shell.sorted_by_attribute('r')
        nshell = len(shell)
        rmin.append( shell[0].r )
        rmax.append( shell[-1].r )
        ravg.append( shell.r.avg() )
        rstd.append( shell.r.std() )
        rmed.append( shell[nshell/2].r )

    return rmin, rmax, ravg, rstd, rmed

def analyse(
        ring,
        planet,
        original_number_of_ringparticles,
        RH          = 1|units.AU,
        boundradii  = False,
        ):
    ring.r = (ring.position - planet.position).lengths()
    ring = ring.sorted_by_attribute('r')

    ring.unbound = np.isnan(orbital_periods(ring, planet).value_in(units.yr))

    escapers    = ring.select(lambda u: u == True, ["unbound"])
    escapers    = escapers.select(lambda r: r > (2|units.AU), ["r"])
    bound       = ring.select(lambda u: u == False, ["unbound"])

    ##FIXME
    #calculate "(surface) density"
    #return radius where density is still higher than???

    bi01 = int(0.01 * len(bound))-1
    bi02 = int(0.02 * len(bound))-1
    bi05 = int(0.05 * len(bound))-1
    bi10 = int(0.10 * len(bound))-1
    bi20 = int(0.20 * len(bound))-1
    bi50 = int(0.50 * len(bound))-1
    bi55 = int(0.55 * len(bound))-1
    bi60 = int(0.60 * len(bound))-1
    bi70 = int(0.70 * len(bound))-1
    bi75 = int(0.75 * len(bound))-1
    bi90 = int(0.90 * len(bound))-1
    bi95 = int(0.95 * len(bound))-1
    bi99 = int(0.99 * len(bound))-1
    bi100 = int(len(bound))-1
    br01 = ring[bi01].r
    br02 = ring[bi02].r
    br05 = ring[bi05].r
    br10 = ring[bi10].r
    br20 = ring[bi20].r
    br50 = ring[bi50].r
    br55 = ring[bi55].r
    br60 = ring[bi60].r
    br70 = ring[bi70].r
    br75 = ring[bi75].r
    br90 = ring[bi90].r
    br95 = ring[bi95].r
    br99 = ring[bi99].r
    br100 = ring[bi100].r

    i01 = int(0.01 * original_number_of_ringparticles)
    i02 = int(0.02 * original_number_of_ringparticles)
    i05 = int(0.05 * original_number_of_ringparticles)
    i10 = int(0.10 * original_number_of_ringparticles)
    i20 = int(0.20 * original_number_of_ringparticles)
    i50 = int(0.50 * original_number_of_ringparticles)
    i55 = int(0.55 * original_number_of_ringparticles)
    i60 = int(0.60 * original_number_of_ringparticles)
    i70 = int(0.70 * original_number_of_ringparticles)
    i75 = int(0.75 * original_number_of_ringparticles)
    i90 = int(0.90 * original_number_of_ringparticles)
    i95 = int(0.95 * original_number_of_ringparticles)
    i99 = int(0.99 * original_number_of_ringparticles)
    r01 = ring[i01].r
    r02 = ring[i02].r
    r05 = ring[i05].r
    r10 = ring[i10].r
    r20 = ring[i20].r
    r50 = ring[i50].r
    r55 = ring[i55].r
    r60 = ring[i60].r
    r70 = ring[i70].r
    r75 = ring[i75].r
    r90 = ring[i90].r
    r95 = ring[i95].r
    r99 = ring[i99].r
    if i50 < len(ring):
        r50 = ring[i50].r
    else:
        r50 = float('inf')
        exit_dissolving()

    outliers = ring.select(lambda r: r > 2*RH, ["r"])
    if boundradii:
        return br10, br20, br50, br55, br60, br70, br75, br90, br95, br99, outliers, escapers
    else:
        return r10, r20, r50, r55, r60, r70, r75, r90, r95, r99, outliers, escapers


def detect_collision(
        ring,
        planet,
        ):
    ring.r = (ring.position - planet.position).lengths()
    ring = ring.sorted_by_attribute('r')
    impacters = ring.select(lambda r: r <= 100*planet.realradius, "r")
    planet.mass += impacters.mass.sum()
    impacters.mass = 0|units.kg

    return impacters



def exit_dissolving():
    logfile.write(
            "#system dissolving\n"
            )
    out_table = open('ring_survival_times.txt','a')
    if thorough:
        mode = "thorough"
    elif quick:
        mode = "quick"
    else:
        mode = "normal"
    out_line = "%s %s %s %s %i %s %s %s %s %s %s %s\n"%(
            planet_mass,
            planet_orbital_period,
            args.fraction_ringsize_hillradius,
            planet_eccentricity,
            pericentre_passages_completed,
            Hillradius.as_quantity_in(units.AU),
            eclipse_time_peri.as_quantity_in(units.day),
            eclipse_time_ap.as_quantity_in(units.day),
            time,
            eclipse_vel_peri.as_quantity_in(units.kms),
            eclipse_vel_ap.as_quantity_in(units.kms),
            mode,
            )
    out_table.write(out_line)
    out_table.close()
    logfile.close()
    errfile.close()
    LRfile.close()
    exit()
    return -1


def find_particle_range(
        particles,
        dmax,
        ):
    try_particle = len(particles)/2
    found = False
    stepsize = len(particles)/4

    while not found:
        if particles[try_particle].r == dmax:
            found           = True
            last_in_shell   = try_particle
        elif particles[try_particle].r > dmax:
            if try_particle == 0:
                found           = True
                last_in_shell   = 0
            elif particles[try_particle-1].r < dmax:
                found           = True
                last_in_shell   = try_particle-1
            else:
                try_particle -= stepsize
                stepsize = max(1, stepsize/2)
        elif particles[try_particle].r < dmax:
            if try_particle == len(particles):
                found           = True
                last_in_shell   = try_particle
            elif particles[try_particle+1].r > dmax:
                found           = True
                last_in_shell   = try_particle
            else:
                try_particle += stepsize
                stepsize = max(1, stepsize/2)
        else:
            found           = True
            last_in_shell   = try_particle
    return last_in_shell + 1

def instantaneous_eccentricity(
        particles,
        planet
        ):

    workingset  = particles.copy()

    m_reduced   = (planet.mass*workingset.mass)/(planet.mass+workingset.mass)

    workingset.v2    = workingset.velocity.lengths_squared()

    #workingset.eps_kinetic   = workingset.specific_kinetic_energy()#0.5 * workingset.v2*workingset.mass
    #workingset.e_potential = workingset.potential()#mass * (-constants.G * planet.mass / workingset.r)

    #workingset.energy      = workingset.mass * (workingset.e_kinetic + workingset.e_potential)
    #workingset.angular_momentum    = workingset.position.cross((workingset.mass.dot(workingset.velocity)))
    
    workingset.epsilon = workingset.specific_kinetic_energy() + workingset.potential()
    workingset.h       = workingset.position.cross(workingset.velocity).lengths()

    workingset.h2      = workingset.h**2

    #print ( 2. * workingset.epsilon * workingset.h2 / ((constants.G*planet.mass)**2) ) 

    mu = (constants.G*planet.mass)**2

    particles.eccentricity    = np.nan_to_num(np.sqrt( 1. + ( 2. * workingset.epsilon * workingset.h2 / mu ) ))
    
    return 


if "__main__" in __name__:
    start = realtime.time()

    snapfile = sys.argv[1]

    allparticles = read_set_from_file(snapfile, 'amuse')

    planet  = allparticles[1]
    
    allparticles.position -= planet.position
    allparticles.velocity -= planet.velocity

    star    = allparticles[0]
    planet  = allparticles[1]
    ring    = allparticles[2:]
    planet_and_ring = allparticles[1:]

    #print star
    #print planet
    
    ring    = ring.sorted_by_attribute('r')

    instantaneous_eccentricity(
            planet_and_ring,
            planet,
            )

    ring.original_r     = ring.r
    ring.original_phi   = ring.phi
    
    #cylindrical_from_xyz(ring)

    ringnr = 0

    nrings = ring.ringnr.max()
    for r in range(nrings):
        p_in_ring = ring.select(lambda x: x == r, ["ringnr"])
        p_in_ring = p_in_ring.sorted_by_attribute('original_phi')
    
        for p in p_in_ring:
            if p.ringnr > ringnr:
                ringnr = p.ringnr
                print("\n")
            try:
                print(p.ringnr, p.original_r.value_in(units.AU), p.r.value_in(units.AU), p.original_phi, p.phi, p.eccentricity)
            except:
                print("-", p.original_r.value_in(units.AU), p.r.value_in(units.AU), p.original_phi, p.phi, p.eccentricity)


    exit()

# coding: utf-8
import logging

from amuse.units import nbody_system, units, constants

from amuse.couple.bridge import CalculateFieldForCodesUsingReinitialize

from amuse.community.hermite0.interface import Hermite
from amuse.community.ph4.interface import ph4
from amuse.community.rebound.interface import Rebound
from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fastkick.interface import FastKick
from amuse.community.seba.interface import SeBa
from amuse.community.kepler_orbiters.interface import Kepler
from amuse.community.bonsai.interface import Bonsai

import numpy as np

logger = logging.getLogger(__name__)

def new_code_gravity(
        converter,
        epsilon,
        p,
        ):
    if p.codes_gravity == "ph4":
        gravity = new_code_gravity_ph4(
                converter,
                epsilon,
                p,
                )
    elif p.codes_gravity == "Hermite":
        gravity = new_code_gravity_hermite(
                converter,
                epsilon,
                p,
                )
    elif p.codes_gravity == "Rebound":
        gravity = new_code_gravity_rebound(
                converter,
                epsilon,
                p,
                )
    else:
        logger.error(
                "No gravity code %s known"%(
                    p.codes_gravity,
                    ),
                )
        exit()
    logger.info(
            "Gravity code %s added"%(
                p.codes_gravity,
                ),
            )
    return gravity
    
def new_code_gravity_rebound(
        converter,
        epsilon,
        p,
        ):
    result  = Rebound(
            converter,
            redirection     = "file",
            redirect_file   = (
                p.dir_codelogs + "/gravity_code.log"
                )
            )
    result.initialize_code()
    print(result.parameters)
    result.parameters.epsilon_squared   = p.particles_epsilon**2
    result.parameters.integrator    = p.codes_gravity_integrator
    result.parameters.timestep      = p.timestep_integrator
    result.parameters.solver        = p.codes_gravity_solver
    result.parameters.boundary      = p.codes_gravity_boundary
    result.parameters.boundary_size = p.codes_gravity_boundary_size
    result.parameters.exact_finish_time = 0
    return result


def new_code_gravity_bonsai(
        converter,
        epsilon,
        p,
        ):
    result  = Bonsai(
            converter,
            redirection     = "file",
            redirect_file   = (
                p.dir_codelogs + "/gravity_code.log"
                )
            )
    result.initialize_code()
    result.parameters.timestep          = p.timestep_integrator
    return result

def new_code_gravity_ph4(
        converter,
        epsilon,
        p,
        ):
    result  = ph4(
            converter,
            redirection     = "file",
            mode            = "gpu" if p.codes_gravity == "gpu" else "cpu",
            redirect_file   = (
                p.dir_codelogs + "/gravity_code.log"
                )
            )
    result.initialize_code()
    result.parameters.epsilon_squared   = epsilon**2
    return result

def new_code_gravity_hermite(
        converter,
        epsilon,
        p,
        ):
    result  = Hermite(
            converter,
            redirection     = "file",
            redirect_file   = (
                p.dir_codelogs + "/gravity_code.log"
                )
            )
    result.parameters.epsilon_squared           = epsilon**2
    result.parameters.end_time_accuracy_factor  = 0
    return result




class kepler_for_bridge(object):
    def __init__(self, converter):
        self.code = Kepler(converter)
        self.model_time = self.code.model_time
        self.central_particle = self.code.central_particle
        self.orbiters = self.code.orbiters
        
    def evolve_model(self, t_end):
        self.code.evolve_model(t_end)
        self.model_time = self.code.model_time

    def get_gravity_at_point(self, radius, x, y, z):
        mass = self.central_particle.mass
        xc, yc, zc = self.central_particle[0].position
        dr2 = ((x-xc)**2+(y-yc)**2+(z-zc)**2+radius**2)
        dr = dr2**0.5
        ax = -mass*(x-xc)/(dr2*dr)
        ay = -mass*(y-yc)/(dr2*dr)
        az = -mass*(z-zc)/(dr2*dr)

        for body in self.orbiters:
            mass = body.mass
            xc, yc, zc = body.position
            dr2 = ((x-xc)**2 + (y-yc)**2 + (z-zc)**2 + radius**2)
            dr = dr2**0.5
            ax -= mass*(x-xc)/(dr2*dr)
            ay -= mass*(y-yc)/(dr2*dr)
            az -= mass*(z-zc)/(dr2*dr)
        ax *= constants.G
        ay *= constants.G
        az *= constants.G
        return ax,ay,az

    def get_potential_at_point(self, radius, x, y, z):
        mass = self.central_particle.mass
        xc, yc, zc = self.central_particle.position
        dr2 = ((x-xc)**2+(y-yc)**2+(z-zc)**2+radius**2)
        dr = dr2**0.5
        phi = -mass/dr

        for body in self.orbiters:
            mass = body.mass
            xc, yc, zc = body.position
            dr2 = ((x-xc)**2+(y-yc)**2+(z-zc)**2+radius**2)
            dr = dr2**0.5
            phi -= mass/dr
        return phi

class non_central_kepler(object):
    #################
    #### Testing ####
    #################
    def __init__(self, converter):
        self.code = Kepler(converter)
        self.model_time = self.code.model_time
        self.central_particle = self.code.central_particle
        self.orbiters = self.code.orbiters
        self.particles = self.orbiters

    #def add_particles(particles):
    #    code.particles.add_particles(particles)
    #
    #def remove_particles():

    def evolve_model(self, t_end):
        self.pc = self.central_particle.position
        self.vc = self.central_particle.velocity

        self.central_particle.position -= self.pc
        self.orbiters.position -= self.pc
        self.central_particle.velocity -= self.vc
        self.orbiters.velocity -= self.vc

        self.code.evolve_model(t_end)
        self.model_time = self.code.model_time

        self.central_particle.position += self.pc
        self.orbiters.position += self.pc
        self.central_particle.velocity += self.vc
        self.orbiters.velocity += self.vc

    def stop(self):
        self.code.stop()

    def get_gravity_at_point(self, radius, x, y, z):
        self.code.get_gravity_at_point(radius, x, y, z)

    def get_potential_at_point(self, radius, x, y, z):
        self.code.get_potential_at_point(radius, x, y, z)


class advance_without_selfgravity(object):
    """
    Code by Lucie Jilkóva
    to advance particles
    """
    def __init__(self, particles, time= 0 |units.Myr):
        self.particles = particles
        self.model_time = time
    
    def evolve_model(self, t_end):
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time= t_end
    
    @property
    def potential_energy(self):
        return quantities.zero
    
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass*self.particles.velocity.lengths()**2).sum()

    def get_gravity_at_point(self, radius, x, y, z):
        fr = (0|units.kg) * constants.G / (1 | units.m**2)
        ax = -0*fr*x/(1|units.m)
        ay = -0*fr*x/(1|units.m)
        az = -0*fr*x/(1|units.m)
        return ax,ay,az


class Star(object):
    def __init__(
            self,
            M = 1.0 | units.MSun,
            R = 1.0 | units.RSun,
            ):
        self.R = R
        self.M = M

    def get_gravity_at_point(self, eps, x, y, z):
        r2  = x**2+y**2+z**2
        r   = r2**0.5
        fr  = constants.G*self.M/r2
        ax  = -fr*x/r
        ay  = -fr*y/r
        az  = -fr*z/r
        return ax,ay,az

    def get_potential_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        c=constant.G*self.M
        phi=c/(r*r2)
        return phi    

    def vcirc(self,r):  
        vc=(constants.G*self.M/r)**0.5
        return vc

def hill_radius(
        eccentricity,
        semimajor_axis,
        minor_mass,
        major_mass,
        ):

    return ((1 - eccentricity) * semimajor_axis * (minor_mass/(3*major_mass))**(1./3) )

def cylindrical_from_xyz(
        particles,
        ):
    particles.r       = particles.position.lengths()
    particles.theta   = np.arccos(particles.z/particles.r)
    particles.phi     = np.arctan2(
            particles.y.value_in(units.AU),
            particles.x.value_in(units.AU),
            )
    return 

def xyz_from_cylindrical(
        orbiters,
        ):
    x = (orbiters.r * np.cos(orbiters.phi))
    y = (orbiters.r * np.sin(orbiters.phi))
    z = 0*x
    return [x,y,z]

def kepler_orbit(
        orbiters, 
        centre, 
        a_over_r = 1.,
        ):
    r = orbiters.position.lengths()
    a = a_over_r * r
    mu = constants.G * centre.mass
    vx +=  (np.sin(orbiters.phi)*(mu * (2./orbiters.r - 1./a))**0.5)
    vy += -(np.cos(orbiters.phi)*(mu * (2./orbiters.r - 1./a))**0.5)
    vz += 0|units.kms
    orbiters.initial_orbital_period = 2 * np.pi * (a**3 / mu)**0.5
    return r, [vx,vy,vz]

def eccentricity(
        orbiters,
        centre,
        ):
    length_unit = units.AU
    speed_unit = units.kms
    rmag = (orbiters.position - centre.position).lengths()
    r = (orbiters.position - centre.position)
    rs = (orbiters.position - centre.position).value_in(length_unit)
    v = (orbiters.velocity - centre.velocity)
    vs = (orbiters.velocity - centre.velocity).value_in(speed_unit)
    h = np.cross(rs,vs)
    mu = (constants.G * centre.mass).value_in(length_unit * speed_unit**2)
    tmp = np.cross(vs,h)/mu
    e = np.array([tmp[:,0]-r[:,0]/rmag,tmp[:,1]-r[:,1]/rmag,tmp[:,2]-r[:,2]/rmag])

    emag = (e[0,:]**2+e[1,:]**2+e[2,:]**2)**0.5

    return e, emag


def orbital_periods(
        allorbiters,
        planet,
        ):
    orbiters = allorbiters
    r   = orbiters.position - planet.position
    v   = orbiters.velocity - planet.velocity
    rabs = r.lengths()
    vabs = v.lengths()
    a   = (constants.G * planet[0].mass) * rabs / (
            (2 * constants.G * planet[0].mass) - rabs*vabs*vabs )
    orbital_periods = (2 * np.pi * (a**3 / (constants.G * planet[0].mass))).sqrt()

    return orbital_periods


def kepler_periastron_velocity(
        mass    = 1.0 | units.MSun,
        epsilon = 0.0,
        a       = 1.0 | units.AU,
        ): 
    mu      = constants.G * mass
    v_per   = np.sqrt( ((1+epsilon)*mu)/((1-epsilon)*a)) 
    return v_per


def kepler_apastron_velocity(
        mass    = 1.0 | units.MSun,
        epsilon = 0.0,
        a       = 1.0 | units.AU,
        ): 
    mu      = constants.G * mass
    v_apo   = np.sqrt( ((1-epsilon)*mu)/((1+epsilon)*a)) 
    return v_apo



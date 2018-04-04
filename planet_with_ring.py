# coding: utf-8
import os,sys

import logging

import numpy as np

import time as clocktime

from amuse.lab import *

from amuse.couple.bridge import Bridge

from amuse.support.console import set_printing_strategy

from parameters import *
from argumentparser import new_simulation_argument_parser
from setup_codes import *
from plotting_new import make_figure
from analysis import run_analysis

def orbital_periods(
        particles,
        centre,
        ):
    r   = particles.position - centre.position
    v   = particles.velocity - centre.velocity
    rabs = r.lengths()
    vabs = v.lengths()
    a   = (constants.G * centre.mass) * rabs / (
            (2 * constants.G * centre.mass) - rabs*vabs*vabs )
    orbital_periods = (2 * np.pi * (a**3 / (constants.G * centre.mass))).sqrt()

    return orbital_periods

class RingSystem(object):

    def __init__(
            self,
            p,
            ):

        self.p  = p
        self.continued_from_backup  = 0

    def setup(self):
        self.setup_codes()
        self.setup_channels()
        
    def initialize_model(self, particles="read"):
        if particles=="read":
            self.particles  = read_set_from_file(p.particles_initial_file, "amuse")
        else:
            self.particles  = particles
        try:
            self.model_time     = self.particles.collection_attributes.model_time
            self.time_offset    = self.particles.collection_attributes.time_offset
            self.p.time_start   = self.model_time
        except:
            self.model_time     = self.p.time_start
            self.time_offset    = 0. | units.yr

        self.define_subgroups()

    def determine_orbit_timestep(self):
        disc_planet_orbital_periods = orbital_periods(
                self.disc,
                self.planets[0],
                )
        try:
            disc_moon_orbital_periods   = orbital_periods(
                    self.disc,
                    self.moons[0],
                    )
        except:
            disc_moon_orbital_periods   = 1|units.Myr
        minimum_disc_planet_period  = np.nanmin(
                disc_planet_orbital_periods.value_in(units.yr),
                ) | units.yr
        minimum_disc_moon_period    = np.nanmin(
                disc_moon_orbital_periods.value_in(units.yr),
                ) | units.yr

        timestep    = min(
                [
                    minimum_disc_planet_period,
                    #minimum_disc_moon_period,
                    ],
                )
        return timestep

    def define_subgroups(self):
        self.stars      = self.particles.select(lambda t: t == "star", ["type"])
        self.planets    = self.particles.select(lambda t: t == "planet", ["type"])
        self.moons      = self.particles.select(lambda t: t == "moon", ["type"])
        self.disc       = self.particles.select(lambda t: t == "disc", ["type"])

    def remove_escapers(self):
        escape_length   = 2.83|units.AU #FIXME: change to Hill radius (for circle orbit?)
        disc        = self.disc.copy()
        disc.d      = (disc.position - self.planets[0].position).lengths()
        escapers    = disc.select(lambda r: r > escape_length, ["d"])

        # Save momentum and energy of escapers
        escaper_momentum    = escapers.total_momentum()
        escaper_kinetic     = escapers.kinetic_energy()
        escaper_potential   = escapers.potential_energy()

        # Remove particles from gravity, and then from the main set
        self.gravity.particles.remove_particles(escapers)
        self.particles.remove_particles(escapers)

        # Re-define subgroups
        self.define_subgroups()

    def setup_codes(self):
        self.gravity    = new_code_gravity(
                p.converter,
                p.particles_epsilon,
                self.p,
                )
        self.gravity.particles.add_particles(self.particles)

        self.collision_detection    = \
                self.gravity.stopping_conditions.collision_detection
        if p.codes_collision_detection:
            self.collision_detection.enable()

        self.gravity.parameters.timestep = (
                p.timestep_accuracy *
                self.determine_orbit_timestep()
                )

        self.old_timestep = self.gravity.parameters.timestep
        print(self.gravity.parameters)

        print(self.gravity.stopping_conditions)

    def setup_channels(self):
        self.from_model_to_gravity  = self.particles.new_channel_to(self.gravity.particles)
        self.from_gravity_to_model  = self.gravity.particles.new_channel_to(self.particles)

    def evolve_model(self, time_end):
        while self.model_time - time_end < - p.timestep_integrator:
            self.time_end   = time_end
            if p.timestep_setting   == "default":
                self.model_time += p.timestep_integrator
            elif p.timestep_setting   == "auto":
                self.new_timestep   = self.gravity.parameters.timestep
                if (
                        round(self.old_timestep.value_in(units.day), 3) != 
                        round(self.new_timestep.value_in(units.day), 3)
                        ):
                    print("New timestep: %s"%self.new_timestep)
                self.old_timestep   = self.new_timestep
                self.model_time += p.timestep
            elif p.timestep_setting   == "fixed":
                self.gravity.parameters.timestep    = p.timestep_integrator
                self.model_time += p.timestep
            elif p.timestep_setting == "orbit":
                self.new_timestep    = (
                        p.timestep_accuracy *
                        self.determine_orbit_timestep()
                        )
                if (
                        round(self.old_timestep.value_in(units.day), 3) != 
                        round(self.new_timestep.value_in(units.day), 3)
                        ):
                    print("New timestep: %s"%self.new_timestep)
                self.gravity.parameters.timestep    = self.new_timestep
                self.old_timestep   = self.new_timestep
                self.model_time += self.gravity.parameters.timestep
            elif p.timestep_setting == "orbit-fixed":
                if self.model_time == 0|units.yr:
                    self.new_timestep    = (
                            p.timestep_accuracy *
                            self.determine_orbit_timestep()
                            )
                    print("New timestep: %s"%self.new_timestep)
                    self.gravity.parameters.timestep    = self.new_timestep
                self.model_time += self.gravity.parameters.timestep
            self.gravity.evolve_model(self.model_time - self.time_offset)
            #print self.gravity.parameters.timestep
            self.model_time = self.gravity.model_time + self.time_offset
            #print self.planets[0].position, self.stars[0].position

            #TODO use stopping conditions to detect/resolve collisions
            if (
                    self.collision_detection.is_set()# and 
                    #(self.gravity.model_time < self.model_time)#FIXME
                    ):
                self.resolve_collision()

        self.sync_model()

        if (
                self.planets[0].position.length().value_in(units.AU) * 0 != 0 or
                self.stars[0].position.length().value_in(units.AU) * 0 != 0
                ):

            print("Something went wrong... Trying to continue from backup %s."%(self.latest_backup))
            self.continue_from_backup()

    def resolve_collision(self):
        print("Resolving collision ", end=' ')
        colliders = Particles()
        colliders.add_particle(self.collision_detection.particles(0))
        colliders.add_particle(self.collision_detection.particles(1))

        # Find out which is the more massive of the two (or pick one if the same)
        for i in range( len( self.collision_detection.particles(0) ) ):
            key0    = self.collision_detection.particles(0)[i].key
            key1    = self.collision_detection.particles(1)[i].key
            p0      = self.particles.select(lambda k: k == key0, ["key"])
            p1      = self.particles.select(lambda k: k == key1, ["key"])
            if p1.mass > p0.mass:
                temp    = p1
                p1      = p0
                p0      = temp
            print("between %s %s and %s %s at t=%s"%(
                    p0.type,
                    p0.key,
                    p1.type,
                    p1.key,
                    self.gravity.model_time,
                    ))

        # Combine masses, calculate    #TODO
        # - new velocities (weighed average)
        # - new positions (weighed average)
        # - new radii (use weighed average of densities of collided particles)

        # In case of massless ring particles: just delete them.
        if p1.mass  == 0 | units.kg:
            self.gravity.particles.remove_particles(p1)
            self.particles.remove_particles(p1)
        else:
            newparticle = Particle()
            newparticle.mass        = colliders.mass.sum()
            newparticle.position    = colliders.center_of_mass()
            newparticle.velocity    = colliders.center_of_mass_velocity()
            #TODO make this better
            newparticle.radius      = p0.radius
            newparticle.type        = p0.type
            newparticle.hoststar    = p0.hoststar
            newparticle.hostplanet  = p0.hostplanet
            newparticle.number      = p0.number

            self.gravity.particles.remove_particle(p0)
            self.gravity.particles.remove_particle(p1)
            self.particles.remove_particle(p0)
            self.particles.remove_particle(p1)

            self.particles.add_particle(newparticle)
            self.gravity.particles.add_particle(newparticle)

        self.define_subgroups()

    def sync_model(self):
        self.from_gravity_to_model.copy_attributes(["x","y","z","vx","vy","vz"])
        self.particles.collection_attributes.time_end   = p.time_end
        self.particles.collection_attributes.model_time = self.model_time
        self.particles.collection_attributes.time_offset = self.gravity.model_time

    def plot(self):
        make_figure(system, p)

    def analysis(self):
        run_analysis(system, p)
        return -1
    
    def save_backup(self):
        self.sync_model()
        backup_name = p.dir_backups+"backup-%syr.hdf5"%(self.model_time.value_in(units.yr))
        write_set_to_file(
                self.particles,
                backup_name,
                'amuse',
                append_to_file  = False,
                )
        self.latest_backup  = backup_name

    def continue_from_backup(self):
        if self.continued_from_backup >= 1: 
            print("Stopping, too many failures")
            exit()
        self.gravity.stop()

        self.particles  = read_set_from_file(
                self.latest_backup,
                'amuse',
                )
        self.stars      = self.particles.select(lambda t: t == "star", ["type"])
        self.planets    = self.particles.select(lambda t: t == "planet", ["type"])
        self.moons      = self.particles.select(lambda t: t == "moon", ["type"])
        self.disc       = self.particles.select(lambda t: t == "disc", ["type"])

        self.model_time     = self.particles.collection_attributes.model_time
        self.time_offset    = self.model_time
        self.setup()

        self.continued_from_backup  += 1



if __name__ in "__main__":
    start_clocktime     = clocktime.time()

    p   = Parameters()
    p   = new_simulation_argument_parser(p)

    particles = read_set_from_file(p.particles_initial_file, "amuse")
    
    p.model_type = particles.collection_attributes.model_type
    p.model_refresh()

    if not os.path.exists(p.dir_simulation):
        os.makedirs(p.dir_simulation)
    if not os.path.exists(p.dir_logs):
        os.makedirs(p.dir_logs)
    if not os.path.exists(p.dir_codelogs):
        os.makedirs(p.dir_codelogs)
    if not os.path.exists(p.dir_plots):
        os.makedirs(p.dir_plots)
    if not os.path.exists(p.dir_backups):
        os.makedirs(p.dir_backups)

    #errfile = open("%s/err.log"%(p.dir_logs),"w")
    #LRfile  = open("%s/LR.log"%(p.dir_logs),"w")
    
    logging.basicConfig(
            level       = logging.DEBUG,
            format      = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt     = "%Y-%m-%d %H:%M:%S %Z",
            filename    = p.dir_logs + "RingSystem.log",
            filemode    = "w",
            )

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    full_formatter   = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    abbr_formatter   = logging.Formatter("%(name)-12s %(levelname)-8s %(message)s")
    simple_formatter = logging.Formatter("%(message)s")
    console.setFormatter(abbr_formatter)

    logger  = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger_filehandler = logging.FileHandler(p.dir_logs + "out.log")
    logger_filehandler.setLevel(logging.DEBUG)
    logger_filehandler.setFormatter(full_formatter)
    logger.addHandler(logger_filehandler)
    logger.addHandler(console)

    analysislogger  = logging.getLogger("analysis")
    analysislogger.setLevel(logging.INFO)
    analysislogger_filehandler  = logging.FileHandler(p.dir_logs + "analysis.log")
    analysislogger_filehandler.setLevel(logging.INFO)
    analysislogger_filehandler.setFormatter(simple_formatter)

    set_printing_strategy(
        "custom", 
        preferred_units = p.log_units_preferred,
        precision       = p.log_units_precision,
        )

    ##### Write descriptions and run information to logfiles
    logging.info("parameters: %s\n"%p)
    #analysislogger.info(
    #        "#time phase pericentre_passages LR10 LR20 LR50 LR55 LR60 LR70 LR75 LR90 LR95 LR99 velocity N_unbound"
    #        )

    system  = RingSystem(p)
    system.initialize_model()
    system.setup()

    time = p.time_start
    
    previous_planet_distance_to_star = (system.planets[0].position - system.stars[0].position).length()
    previous_planet_orbital_phase = np.arctan2(
            (system.planets[0].y - system.stars[0].y).value_in(units.AU),
            (system.planets[0].x - system.stars[0].x).value_in(units.AU),
            )

    drift_position = 0 * system.stars[0].position
    drift_velocity = 0 * system.stars[0].velocity

    cleanup_step                            = 0
    figure_num                              = 0
    backup_num                              = 0
    previous_pericentre_passages_completed  = 0
    pericentre_passages_completed           = 0
    apocentre_passages_completed            = 0

    next_plot_time                          = time
    next_analyse_time                       = time
    next_backup_time                        = time

    system.sync_model()

    system.save_backup()

    planet_orbital_phase = np.arctan2(
        (system.planets[0].y - system.stars[0].y).value_in(units.AU),
        (system.planets[0].x - system.stars[0].x).value_in(units.AU),
        )

    print(time, system.model_time, system.gravity.model_time, planet_orbital_phase, (system.planets[0].position-system.stars[0].position).length(), len(system.particles), len(system.gravity.particles), (system.planets[0].velocity-system.stars[0].velocity).length(), system.disc.y.max()- system.disc.y.min(), (system.disc.y.max()-system.disc.y.min())/(system.planets[0].vy-system.stars[0].vy))

    system.plot()
    ##################################################
    ######   Main time loop starts here         ######
    ##################################################
     
    while time < p.time_end:
        time += p.timestep

        system.evolve_model(time)
        system.sync_model()

        planet_distance_to_star = (system.planets[0].position - system.stars[0].position).length()
        planet_orbital_phase = np.arctan2(
                (system.planets[0].y - system.stars[0].y).value_in(units.AU),
                (system.planets[0].x - system.stars[0].x).value_in(units.AU),
                )

        previous_planet_distance_to_star    = planet_distance_to_star
        previous_planet_orbital_phase       = planet_orbital_phase
        
        center          = system.planets[0].position - system.stars[0].position
        center_velocity = system.planets[0].velocity - system.stars[0].velocity
        dist            = center.length().as_quantity_in(units.AU)
        vel             = center_velocity.length().as_quantity_in(units.kms)

        if time > next_plot_time:
            next_plot_time += p.timestep_plotting
            system.plot()

        if time > next_analyse_time:
            next_analyse_time += p.timestep_analysis
            system.analysis()
        
        if time > next_backup_time:
            next_backup_time += p.timestep_backup
            system.save_backup()

        #print time, system.model_time, system.gravity.model_time, planet_orbital_phase, (system.planets[0].position-system.stars[0].position).length(), len(system.particles), len(system.gravity.particles), (system.planets[0].velocity-system.stars[0].velocity).length(), system.disc.y.min(), system.disc.y.max()
        
        #system.remove_outliers()

        system.remove_escapers()

        print(time, system.model_time, system.gravity.model_time, planet_orbital_phase, (system.planets[0].position-system.stars[0].position).length(), len(system.particles), len(system.gravity.particles), (system.planets[0].velocity-system.stars[0].velocity).length(), system.disc.y.max()- system.disc.y.min(), (system.disc.y.max()-system.disc.y.min())/(system.planets[0].vy-system.stars[0].vy))
            

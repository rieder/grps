import numpy as np

from amuse.units import units, nbody_system, constants
from amuse.support.console import set_printing_strategy
from setup_codes import hill_radius

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

class Parameters(object):

    def __init__(self, model_type="BR80"):
        ### General #####################

        ### Initial conditions ##########
        self.model_type                 = model_type

        ### Codes and code settings #####
        #self.codes_gravity              = "ph4"
        self.codes_gravity              = "Rebound"
        self.codes_gravity_gpu          = False
        #self.codes_gravity_integrator   = "ias15"
        #self.codes_gravity_integrator   = "hybrid"
        self.codes_gravity_integrator   = "whfast"
        #self.codes_gravity_integrator   = "wh"
        #self.codes_gravity_integrator   = "sei"
        #self.codes_gravity_integrator   = "leapfrog"
        #self.codes_gravity_solver       = "tree"
        self.codes_gravity_solver       = "compensated"
        self.codes_gravity_boundary     = "open"
        self.codes_gravity_boundary_size= 500000 | units.AU
        #self.codes_bridge_feedback      = "FastKick"
        #self.codes_bridge_feedback_gpu  = False
        #self.codes_sph                  = "gadget"
        #self.codes_sph_nworkers         = 4
        #self.codes_collision_detection  = True
        self.codes_collision_detection  = False

        ### Particles ###################
        self.particles_initial_file     = "InitialConditions.hdf5"
        #self.particles_model            = "uniform" # or random

        ### Timesteps ###################
        self.timestep                   = 5.5  | units.yr
        self.timestep_analysis          = 5.5  | units.yr  
        self.timestep_plotting          = 0.1  | units.yr  
        self.timestep_backup            = 60.5 | units.yr

        self.timestep_setting           = "auto"
        #self.timestep_setting           = "default"
        #self.timestep_setting           = "fixed"
        #self.timestep_setting           = "orbit"
        #self.timestep_setting           = "orbit-fixed"
        ## Number of steps per orbit/nbody time unit (should be 0.01 or so)
        self.timestep_accuracy          = 0.01
        #self.timestep_integrator        = 0.002 | units.yr
        
        ### Simulation duration #########
        self.time_start                 =   0.0 | units.yr
        self.time_end                   = 1.0e5 | units.yr

        ### Plotting ####################
        #self.plot_bins                  = 128
        self.plot_dpi                   = 150
        self.plot_figsize               = (3,3)
        self.plot_minx                  = -2.0 | units.AU
        self.plot_maxx                  =  2.0 | units.AU
        self.plot_miny                  = -2.0 | units.AU
        self.plot_maxy                  =  2.0 | units.AU
        self.plot_bgcolor               = "white"
        self.plot_axes_x                = "x"
        self.plot_axes_y                = "y"
        self.plot_tics                  = True
        self.plot_fontsize              = 20
        #self.plot_colormap              = ""
        #self.plot_vmin                  = 1.2
        #self.plot_vmax                  = 4.0

        ### Units #######################
        self.units_time                 = units.yr
        self.units_length               = units.km
        self.units_speed                = units.kms
        self.units_mass                 = units.MJupiter

        ### Logging #####################
        self.log_units_preferred        = [
                units.MJupiter, 
                units.AU, 
                units.yr, 
                units.kms,
                ]
        self.log_units_precision        = 5
        set_printing_strategy(
                "custom", 
                preferred_units         = self.log_units_preferred, 
                precision               = self.log_units_precision,
                )

        self.model_refresh()
        self.variable_refresh()

    def model_refresh(self):
        # Here go parameters that (may) depend on the model_type

        ### Initial conditions ##########
        self.seed                       = 42
        self.stars                      = 1
        self.planets                    = 1
        self.moons                      = 0
        self.discs                      = 1

        self.stars_mass                 = [
                1.0 | units.MSun,
                ]
        self.stars_radius               = [
                1.0 | units.RSun,
                ]
        self.stars_position             = [
                np.array([0,0,0])|units.AU,
                ]
        self.stars_velocity             = [
                np.array([0,0,0])|units.kms,
                ]
        self.stars_eccentricity         = [
                0.0,
                ]
        self.stars_semimajoraxis        = [
                0.0 | units.AU,
                ]
        self.planets_host               = [
                0,
                ]
        self.planets_mass               = [
                float(self.model_type[2:]) | units.MJupiter,
                #80.0 | units.MJupiter,
                ]
        self.planets_radius             = [
                4.0 | units.RJupiter,
                ]
        self.planets_eccentricity       = [
                0.70 if self.model_type[0]=="A" else 
                0.65 if self.model_type[0]=="B" else 
                0.60 if self.model_type[0]=="C" else 
                0.0,
                ]
        self.planets_period             = [
                11 | units.yr,
                ]
        self.planets_semimajoraxis      = [
                semimajoraxis(self.stars_mass[0], self.planets_mass[0], self.planets_period[0]),#5.06921659724 | units.AU,
                ]
        if self.moons == 1:
            self.moons_host                 = [
                    0,
                    ]
            self.moons_mass                 = [
                    6.417e23 | units.kg,
                    ]
            self.moons_radius                = [
                    0.25 | units.RJupiter,
                    ]
            self.moons_eccentricity         = [
                    0.00,
                    ]
            self.moons_semimajoraxis        = [
                    0.2 | units.AU,
                    ]
            self.moons_direction            = [
                    "retrograde" if self.model_type[1]=="R" else
                    "prograde",
                    ]
        else:
            self.moons_host                 = [
                    ]
            self.moons_mass                 = [
                    ]
            self.moons_radius                = [
                    ]
            self.moons_eccentricity         = [
                    ]
            self.moons_semimajoraxis        = [
                    ]
        self.discs_host                 = [
                ["planet",0],
                ]
        self.discs_n                    = [
                5000,
                ]
        self.discs_rmax                 = [
                hill_radius(
                    self.planets_eccentricity[0],
                    self.planets_semimajoraxis[0],
                    self.planets_mass[0],
                    self.stars_mass[0],
                    )#0.50 | units.AU,
                ]
        self.discs_rmin                 = [
                0.25 * self.discs_rmax[0],
                ]
        self.discs_density              = [
                3 | units.g*units.cm**-3,
                ]
        self.discs_mass                 = [
                0 | units.MJupiter,
                ]
        self.discs_zrange               = [
                1.0 | units.RJupiter,
                ]
        self.discs_model                = [
                #"random",
                "powerlaw",
                ]
        self.discs_direction            = [
                "retrograde" if self.model_type[1]=="R" else
                "prograde",
                ]
        self.discs_nrings               = [
                50,
                ]# in case of non-random disk


        ### Converter ###################
        self.scale_M                    = self.planets_mass[0]#80.0  | units.MJupiter
        self.scale_R                    = self.discs_rmin[0]#1.00  | units.AU
        self.converter                  = nbody_system.nbody_to_si(
                self.scale_M,
                self.scale_R,
                )

        ### Timesteps ###################
        self.timestep_integrator        = \
                self.converter.to_si(self.timestep_accuracy | nbody_system.time)

        ### Directories #################
        self.model_name                  = "%s%s-%s"%(
                self.model_type[0],
                self.model_type[2:],
                self.discs_direction[0],
                )

        self.dir_simulation             = "./%s/"%self.model_name
        self.dir_logs                   = self.dir_simulation + "logs/"
        self.dir_codelogs               = self.dir_simulation + "codelogs/"
        self.dir_plots                  = self.dir_simulation + "plots/"
        self.dir_backups                = self.dir_simulation + "backups/"


        ### Particles ###################
        self.particles_epsilon          = 0. | units.RJupiter 
        ## When interaction moon-disc is important, this should be small:
        ## not larger than the size of the smallest moon
        ## Ideally, this would be set on a per-particle basis... Or need to handle collisions.
        #self.particles_interaction_radius   = 100 | units.AU


    def variable_refresh(self):
        # (Re-)initialise parameters that depend on another parameter
        # This function should be called after any of these parameters
        # changes!

        ### Directories #################
        self.dir_initialconditions      = self.dir_simulation + "ICs/"
        self.dir_logs                   = self.dir_simulation + "logs/"
        self.dir_codelogs               = self.dir_simulation + "code_logs/"
        self.dir_plots                  = self.dir_simulation + "plots/"

        ### Logging #####################
        self.log_output                 = self.dir_logs + "output.log"
        self.log_energy                 = self.dir_logs + "energy.log"
        #self.log_encounters             = self.dir_logs + "encounters.log"


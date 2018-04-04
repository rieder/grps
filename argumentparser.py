import argparse

def new_simulation_argument_parser(p):

    parser = argparse.ArgumentParser()

    parser.add_argument(
            '-g',
            dest    = 'codes_gravity',
            type    = str,
            default = p.codes_gravity,
            help    = "Gravity code to use [%s]"%p.codes_gravity,
            )
    parser.add_argument(
            '-a',
            dest    = 'timestep_analysis',
            type    = float,
            default = p.timestep_analysis.value_in(p.units_time),
            help    = "Analysis timestep [%s] (0: disabled)"%p.timestep_analysis,
            )
    parser.add_argument(
            '-p',
            dest    = 'timestep_plotting',
            type    = float,
            default = p.timestep_plotting.value_in(p.units_time),
            help    = "Plotting timestep [%s] (0: disabled)"%p.timestep_plotting,
            )
    parser.add_argument(
            '-b',
            dest    = 'timestep_backup',
            type    = float,
            default = p.timestep_backup.value_in(p.units_time),
            help    = "Backup timestep [%s] (0: disabled)"%p.timestep_backup,
            )
    parser.add_argument(
            '-t',
            dest    = 'timestep',
            type    = float,
            default = p.timestep.value_in(p.units_time),
            help    = "Integration timestep [%s] (0: auto)"%p.timestep_backup,
            )
    parser.add_argument(
            '-e',
            dest    = 'time_end',
            type    = float,
            default = p.time_end.value_in(p.units_time),
            help    = "End time [%s]"%p.time_end,
            )
    parser.add_argument(
            '-i',
            dest    = 'particles_initial_file',
            type    = str,
            default = p.particles_initial_file,
            help    = "Use this model [%s]"%p.particles_initial_file,
            )

    args = parser.parse_args()

    p.codes_gravity             = args.codes_gravity
    p.timestep_analysis         = args.timestep_analysis | p.units_time
    p.timestep_plotting         = args.timestep_plotting | p.units_time
    p.timestep_backup           = args.timestep_backup | p.units_time
    p.timestep                  = args.timestep | p.units_time
    p.time_end                  = args.time_end | p.units_time
    p.particles_initial_file    = args.particles_initial_file

    return p

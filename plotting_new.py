# coding: utf-8
import os,sys
import argparse
import numpy as np
import time as clocktime

from amuse.lab import *

from amuse.community.kepler_orbiters.interface import Kepler
from amuse.couple.bridge import Bridge
from amuse.support.console import set_printing_strategy

from parameters import Parameters

def new_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-l',
            dest    = 'savefile',
            type    = str,
            default = "none",
            help    = "Load this savefile",
            )

    args = parser.parse_args()
    return args


def make_figure(
        system,
        p,
        center          = np.array([0, 0, 0])|units.AU,
        figfile         = "plot-%011.5fyr.png",
        dpi             = 150,
        alpha           = 0.75,
        fontsize        = 24,
        tight           = False,
        plot_spines     = True,
        darkmode        = True,
        minimalmode     = True,
        ):

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    if darkmode:
        plt.style.use('dark_background')


    stars   = system.stars
    planets = system.planets
    moons   = system.moons
    disc    = system.disc

    center = planets[0].position - stars[0].position
    
    pointsize = 2.0

    fig         = plt.figure(figsize=(5,5),dpi=dpi)
    ax          = fig.add_subplot(111, aspect='equal')

    if minimalmode:
        plot_spines = False

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        tight = True
    else:
        ax.set_xlabel(p.plot_axes_x, fontsize=p.plot_fontsize)
        ax.set_ylabel(p.plot_axes_y, fontsize=p.plot_fontsize)

    ax.spines["top"].set_visible(plot_spines)
    ax.spines["right"].set_visible(plot_spines)
    ax.spines["bottom"].set_visible(plot_spines)
    ax.spines["left"].set_visible(plot_spines)
    

    boundaries  = [ 
            center[0].value_in(p.log_units_preferred[1])+p.plot_minx.value_in(p.log_units_preferred[1]), 
            center[0].value_in(p.log_units_preferred[1])+p.plot_maxx.value_in(p.log_units_preferred[1]), 
            center[1].value_in(p.log_units_preferred[1])+p.plot_miny.value_in(p.log_units_preferred[1]), 
            center[1].value_in(p.log_units_preferred[1])+p.plot_maxy.value_in(p.log_units_preferred[1]),
            ]

    if darkmode:
        particlecolor = "white"
        ax.set_facecolor("black")
    else:
        particlecolor = "black"
        ax.set_facecolor("white")

    ax.scatter(
            (disc.x - stars[0].x).value_in(p.log_units_preferred[1]),
            (disc.y - stars[0].y).value_in(p.log_units_preferred[1]),
            s           = pointsize,
            color       = particlecolor,
            alpha       = alpha,
            marker      = "o",
            lw          = 0,
            edgecolors  = "none",
            )

    ax.axis( boundaries )

    starcircles = []
    for star in stars:
        starcircles.append(plt.Circle(
                (
                    (star.x - stars[0].x).value_in(p.log_units_preferred[1]), 
                    (star.y - stars[0].y).value_in(p.log_units_preferred[1]),
                    ),
                10* star.radius.value_in(p.log_units_preferred[1]),
                color = 'yellow',
                ))
    for starcircle in starcircles:
        fig.gca().add_artist(starcircle)

    planetcircles = []
    for planet in planets:
        planetcircles.append(plt.Circle(
                (
                    (planet.x - stars[0].x).value_in(p.log_units_preferred[1]), 
                    (planet.y - stars[0].y).value_in(p.log_units_preferred[1]),
                    ),
                10* planet.radius.value_in(p.log_units_preferred[1]),
                facecolor   = 'orange',
                edgecolor   = "none",
                fill        = True,
                ))
    for planetcircle in planetcircles:
        fig.gca().add_artist(planetcircle)

    mooncircles = []
    for moon in moons:
        mooncircles.append(plt.Circle(
                (
                    (moon.x - stars[0].x).value_in(p.log_units_preferred[1]), 
                    (moon.y - stars[0].y).value_in(p.log_units_preferred[1]),
                    ),
                10*moon.radius.value_in(p.log_units_preferred[1]),
                facecolor   = 'red',
                edgecolor   = "none",
                fill        = True,
                ))
    for mooncircle in mooncircles:
        fig.gca().add_artist(mooncircle)

    fig.savefig(
            p.dir_plots + figfile%(system.model_time.value_in(units.yr)),
            frameon     = False,
            dpi         = dpi,
            bbox_inches = 'tight' if tight else None,
            )
    plt.close(fig)

    return

if __name__ in "__main__":

    p = Parameters()
    start_clocktime     = clocktime.time()
    args                = new_argument_parser()
    verbose             = False

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    #from colorschemes_pault import *

    particles = read_set_from_file(args.savefile, 'amuse')
    particles.position -= particles[1].position
    star = particles[0]
    planet = particles[1]
    moons = particles[2:3]
    #moons = Particles()
    ring = particles[2:]

    center = planet.position

    make_figure(
            ring,
            0,
            star,
            planet,
            moons,
            bins            = 512,
            center          = center,
            figfile         = args.savefile+"%s-%i.png",
            figname         = "standard",
            time            = False,
            method          = "scatter",
            plot_hillradius = False,
            plot_gravity    = False,
            color           = "black",
            bgcolor         = "white",
            alpha           = 1.0,#0.75,
            plot_spines     = False,
            xlabel          = "X [%s]"%(p.log_units_preferred[1]),
            ylabel          = "Y [%s]"%(p.log_units_preferred[1]),
            pointsize       = 1,
            tight           = True,
            #linewidth       = 4,

            #hillradius      = (
            #    ( (planet.x-star.x)**2+(planet.y-star.y)**2 )**0.5 * \
            #            (planet.mass/(3*star.mass))**(1./3) ),
            #drawtime        = True,
            )
    #make_figure(
    #        ring,
    #        figure_num,
    #        star,
    #        planet,
    #        moons,
    #        bins            = 1440,
    #        center          = star[0].position,
    #        figfile         = figdir+"zoomout/"+figfile,
    #        figname         = "zoomout",
    #        time            = time,
    #        method          = "scatter",
    #        plot_gravity    = True,
    #        convert_nbody   = ring_converter,
    #        )
    #make_figure(
    #        ring,
    #        figure_num,
    #        star,
    #        planet,
    #        moons,
    #        bins    = 512,
    #        center  = center,
    #        figfile = figdir+"zoomin/"+figfile,
    #        figname = "zoomin",
    #        time    = time,
    #        method  = "scatter"
    #        )

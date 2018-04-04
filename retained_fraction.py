import os,sys
from amuse.lab import *
import numpy as np

initial_particles = read_set_from_file(sys.argv[1], "amuse")
final_particles   = read_set_from_file(sys.argv[2], "amuse")

initial_disc      = initial_particles.select(lambda x: x == "disc", ["type"])
initial_planet    = initial_particles.select(lambda x: x == "planet", ["type"])[0]
initial_star      = initial_particles.select(lambda x: x == "star", ["type"])[0]
initial_disc.position -= initial_planet.position
initial_disc.velocity -= initial_planet.velocity
initial_planet.position -= initial_star.position
initial_planet.velocity -= initial_star.velocity

final_disc        = final_particles.select(lambda x: x == "disc", ["type"])
final_planet      = final_particles.select(lambda x: x == "planet", ["type"])[0]
final_star        = final_particles.select(lambda x: x == "star", ["type"])[0]
final_disc.position -= final_planet.position
final_disc.velocity -= final_planet.velocity
final_planet.position -= final_star.position
final_planet.velocity -= final_star.velocity

final_disc.distance     = final_disc.position.lengths()
initial_disc.distance   = initial_disc.position.lengths()
final_disc              = final_disc.sorted_by_attribute("distance")
initial_disc            = initial_disc.sorted_by_attribute("distance")

initial_disc_set  = set(initial_disc.key)
final_disc_set    = set(final_disc.key)

retained_fraction = initial_particles.select(lambda x: x in final_disc_set, ["key"])

#retained_fraction.position -= initial_planet.position

hill              = initial_disc.position.lengths().max()
most_distant_init = retained_fraction.distance.max()
most_distant_final= final_disc.distance.max()
x = int(0.99*len(retained_fraction))
r99_distant_init  = retained_fraction[x].distance
r99_distant_final = final_disc[x].distance

final_disc              = final_disc.sorted_by_attribute("y")
xx = int(0.005 * len(retained_fraction))
dy = final_disc[-1].y - final_disc[0].y
dy99 = final_disc[-1-xx].y - final_disc[xx].y

print(sys.argv[1], end=' ') 
print(hill.as_quantity_in(units.AU), end=' ')
print(most_distant_init/hill, most_distant_init.as_quantity_in(units.AU), end=' ')
print(r99_distant_init/hill,  r99_distant_init.as_quantity_in(units.AU), end=' ')
print(most_distant_final/hill, most_distant_final.as_quantity_in(units.AU), end=' ')
print(r99_distant_final/hill, r99_distant_final.as_quantity_in(units.AU), end=' ')
print(dy/hill, dy.as_quantity_in(units.AU), end=' ')
print(dy99/hill, dy99.as_quantity_in(units.AU), end=' ')
print(final_planet.x.as_quantity_in(units.AU), final_planet.y.as_quantity_in(units.AU), end=' ')
print(final_planet.vx.as_quantity_in(units.kms), final_planet.vy.as_quantity_in(units.kms), end=' ')
print((dy/final_planet.vy).as_quantity_in(units.day), (dy99/final_planet.vy).as_quantity_in(units.day))

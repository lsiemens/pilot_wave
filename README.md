# Pilot Wave
Visualize 3d quantum wave functions using pilot wave theory and the Madelung
equations. This visualization is based on the video *A Better Way To Picture
Atoms* by MinutePhysics. In 2D a spatial quantum wave function can be fully
represented, the wave function is described by its amplitude and phase which can
be represented at each point by HSL color (or an equivalent system) where the
amplitude is represented by the luminance and the phase is represented by the
hue at each point. In 3D wave functions are often represented by a level
surface, sometimes colored using hue to represent the phase. These
representations do express or attempt to indicate the form of the phase and
amplitude, but this does not fully capture the volumetric structure of the
orbitals and how the phase relates to the orbital angular momentum of the state.

In the same way that fog illuminated by a laser can show the structure and
motion of air currents, test particles scattered through the domain can be used
to give a different representation of a wave function that captures some aspects
of the volumetric structure and can be related to the orbital angular momentum.
The core of the visualization is to distribute test particles such that the
number density of test particles is proportional to the squared amplitude of the
wave function, and give velocities to the test particles to match the
probability current of the wave function. Then the probability of observing the
quantum particle at a point is show by the density of test particles and the
orbital angular momentum is represented by orbital motion of the test particles.
Given standard interpretations of quantum mechanics the test particles are
incompatible with the properties of quantum objects (they have definite position
and momentum/velocity) they simply are an abstract representation of the
probability of measuring the state to be at a point and provide a representation
of angular momentum. However, given the De Broglie-Bohm formulation (pilot wave
theory) these test particles and their tracks represent possible physical paths
of particles, in this case the motion of the test particles is not just an
analogy for more abstract properties but is real physical motion.

The mathematics is based on the Madelung equations which are related to the
equations of pilot wave theory and produce identical particle tracks. 

Represent a collection of particles but possible locations where the quantum
particle could be observed. Also the motion of the test particles 

This code is based on the [tutorial](http://www.opengl-tutorial.org/beginners-tutorials/tutorial-1-opening-a-window/)
from https://github.com/opengl-tutorials/ogl.

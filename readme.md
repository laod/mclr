This is a coffeescript implementation of Florian Boesch's [Minecraft Like Rendering Experiments](http://codeflow.org/entries/2010/dec/09/minecraft-like-rendering-experiments-in-opengl-4/). Thanks to him for the inspiration and code.

Warning
-------

This code has only been tested on a recent chrome in linux. It's also VERY SLOW right now. I have a long way to go. Open your js console and be patient.

Background
----------

The idea was to learn about webgl, glsl, and gpgpu.

Interesting diversions from the original:

 * no geometry shaders, sadly
 * webgl gpgpu (ambient occlusion is calculated in shaders and read back, as opposed to being pure shader)

Makes use of three.js, underscore.js, and jquery.

Requirements
------------

 * Browser with webgl with fp texture extension support

To-do
-----

 * Refactor and clean up
 * Texturing
 * Observer light
 * Improve ao parallelism and move to progressive refinement
 * Gamma correction?
 * understand windows js memory limits (s=128, ts=2048)
 * understand why shaders break at s=128

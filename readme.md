This is a coffeescript implementation of Florian Boesch's [Minecraft Like Rendering Experiments](http://codeflow.org/entries/2010/dec/09/minecraft-like-rendering-experiments-in-opengl-4/). Thanks to him for the inspiration and code.

Warning
-------

This code has only been tested on a recent chrome in linux. It's also VERY SLOW right now. I have a long way to go. Open your js console and be patient.

Background
----------

The idea was to learn about webgl and glsl.

Interesting diversions from the original:

 * no geometry shaders, sadly
 * webgl gpgpu (ambient occlusion is calculated in shaders)

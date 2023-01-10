---
layout: page
title: contours
description: clever visualization of 2D Perlin noise
img: assets/img/contour_preview.png
importance: 1
category: fun
---

As part of my new years resolutions, I am aiming to become proficient in a programming language entirely removed from my primary PhD research. Part of this is also to endeavour into a world that allows me to be more creative, let's say. This project is my first attempt to do so.

Perlin noise is a well-utilized algorithm to procedurally generate gradient noise such that the resulting $n$-dimensional noise landscape is quasi-random, but with constraints such that it appears continuous and thus life-like. Originally developed to procedurally build the world of *Tron* (1982), this noise is the *de facto* standard for video games and CGI landscape development. Since there is so much knowledge upon which we can further develop, this seems an appropriate place to start in my creative ventures. 

If I'm going to learn a new programming language, I'm going to go all in. I aimed to create a piece of generative art, such that a string can be used to seed a random number generator, and so repeating the generation with the same string produces the same piece of art. This is easily done by first hashing the seed string into a 128-bit hash which then seeds the *simple fast counter*, for which plenty of documentation exists. From here, we make the art!

Generating 2D perlin noise generates an image $$\mathcal{I}(\mathsf{x},\mathsf{y})\in [0,1)\times[0,1)\subset\mathbb{R}^2\;\forall (\mathsf{x},\mathsf{y})\subset[\mathsf{x}_\mathsf{min},\mathsf{x}_\mathsf{max}]\times[\mathsf{y}_\mathsf{min},\mathsf{y}_\mathsf{max}]\sim\mathcal{R}$$, where the bounds and spacing of points inside the rectangle are determined by the random number generator. In order to create a visually appealing image, we round the images values to discrete steps to generate clear *level sets*. Here, a level set is the region of the image of constant value between the contours:

$$
\eta_C\equiv\left\{(\mathsf{x},\mathsf{y})\in\mathcal{R} \;|\;\mathcal{I}(\mathsf{x},\mathsf{y})=C\right\}
$$

From here, we iterate over the contours, inside each we create sets of randomly generated points such that the density of points is uniform over $$\mathcal{I}$$. After asssigning colormaps (also randomly assigned per seed string), we get images like the one given below.

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% responsive_image path: assets/img/contour_example.png title: "example image" class: "img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
    </div>
</div>
<div class="caption">
   An example of the patterns generated from the algorithm. This image uses the seed <code class="language-plaintext highlighter-rouge">tristan</code>
</div>

If you want to try different seeds to generate your own patterns, feel free to try it out [here](../contours.html). By default, the image rendered will use my name. To pass different seeds, pass `?seed=<characters_you_want_to_use>` to the URL. If you want a higher resolution render (for a poster for example), you can pass a second argument with the `&` symbol, and use `density=<how_dense_you_want_the_pixels>`, where you can use densities $$<$$ 4. For example, to use the word `January` with as dense a render as possible, add `?seed=January&density=4` to the URL. 

Happy rendering :)
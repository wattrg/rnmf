# rnmf
Numerical method framework written in Rust.


## The idea
My goal is to implement these eventually, and not necessarily in this order
* Support many types of mesh in 1, 2, and 3D
   * Cartesian (eventually with ability to adaptively refine)
   * Structured, non-cartesian
   * Unstructured
* Provide vectorised/parallelisable data structures for the types of mesh (e.g. tile-able data structures)
* Provide framework for simple PDEs with FD discretisation
* Provide support for particles
* Use rLua to be able to use Lua to configure the simulation. Ideally will be able to run a sim in two ways:
    * Passing Lua script to executable as command line argument
    * Running the simulation directly from lua script

![rotator_and_cutaway.jpg](https://bitbucket.org/repo/Rp975y/images/1513682443-rotator_and_cutaway.jpg)

# What is Rayleigh? #

Rayleigh is a 3-D convection code designed for the study of planetary and stellar dynamos.  In particular, it evolves the incompressible and anelastic MHD equations in spherical geometry using a pseudo-spectral approach.  Rayleigh employs spherical harmonics in the horizontal direction and Chebyshev polynomials in the radial direction.  Rayleigh has been developed with NSF support through the Computational Infrastructure for Geodynamics ([CIG](https://geodynamics.org/cig/news/newsletters/may-2016/)).

# What Makes Rayleigh Unique? #
The pseudo-spectral nature of Rayleigh means that its parallelization necessarily relies heavily on **global communication** patterns.  That said, Rayleigh's parallelization is based around a 2-D domain decomposition and large-message-size all-to-alls.   The end result is a pseudo-spectral code optimized for petascale machines.  Rayleigh has demonstrated highly efficient strong scaling on  131,000 cores of the Mira Blue Gene/Q supercomputer for problems with approximately 2048^3 grid point (2048 spherical harmonics).  Performance numbers from Mira are shown below.
 
![mira_performance.jpg](https://bitbucket.org/repo/Rp975y/images/3897197863-mira_performance.jpg)


# When Will Rayleigh Be Released?#
Rayleigh is sleighted for release on June 19, 2016.  In the meantime, please pardon the mess.  The code and (and presently incomplete) documentation are in a state of flux.  See the Wiki and the documentation directory (in the source) for tips.  Some videos corresponding to the images below can be found [here](http://www.youtube.com/user/feathern24).
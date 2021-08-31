## DESCRIPTION

galbispectra is a code built on top of the SONG (Second Order Non-Gaussianity) - see https://github.com/coccoinomane/song.git 

SONG is a second-order Boltzmann code which computes the non-linear evolution of the Universe in order to predict cosmological observables such as the bispectrum of the Cosmic Microwave Background (CMB). The module galbispectra is an extension to SONG which allows for the computation of the angular galaxy bispectrum.

More precisely, SONG (Second Order Non-Gaussianity) is a second-order Boltzmann code, as it solves the Einstein and Boltzmann equations up to second order in the cosmological perturbations; the physics, mathematics and numerics of SONG are described extensively in G. Pettinari's publicly available [PhD thesis][10], especially in Chapters 5 and 6.

The main additions with respect to SONG are modifications or the addition of the following files

- galbispectra2.c
- galbispectra2.h
- song.c
- input2.c
- input2.h

## GETTING STARTED
1. First download the and install the WIGXJPF library.

http://fy.chalmers.se/subatom/wigxjpf/

2. Download galbispectra. You can do so either from the project's page <https://github.com/sedlawrence/galbispectra/> or from the command line with:

    git clone --recursive https://github.com/sedlawrence/galbispectra.git
    
    The `--recursive` flag is important because it allows to download both SONG and CLASS in one command.

3. Edit the makefile such that it reads (i.e. ensure that the wigxjpf directories are linked to the locations on your system):

	# Header files and libraries
	INCLUDES = -I../include -I../$(CLASS_DIR)/include -I/<YOUR-PATH-TO>/wigxjpf-1.11/inc/
	LDFLAGS = -lm -lwigxjpf

	CFLAGS += -I/<YOUR-PATH-TO>/wigxjpf-1.11/inc/
	LDFLAGS += -L/<YOUR-PATH-TO>/wigxjpf-1.11/lib/
	
4. To compile and make a test run of SONG, enter SONG's directory and execute the following commands from your terminal:

    make song
    ./song explanatory_template.ini pre/quick_song_run.pre
	
    I have been working with the compiler gcc --version gcc (GCC) 10.2.1 20200804 (Red Hat 10.2.1-2)

## MAC SUPPORT

On Mac Os X the default `gcc` compiler is `clang`, which does not natively support OpenMP.
To compile SONG without parallel support, comment the lines with `-fopenmp` in the Makefile.
If instead you want parallel support, I suggest either one of the following:

* Use the GNU `gcc` compiler, which can be downloaded from [Macports], [Homebrew] or [HPC].
* Install the clang-optimized version of OpenMP (see for ex. https://stackoverflow.com/a/39843038/2972183).


## DIRECTORY STRUCTURE
The directory structure of SONG is important to learn how the code works:

* The 'source' directory contains the main source files in C. Each file corresponds to a module in SONG.

* The 'tools' directory contains accessory source files in C with purely numerical functions or utility functions.

* The 'main' directory contains the main source files, i.e. the executable files, including song.c.

* The 'python' directory contains Python scripts to inspect some of SONG outputs (credits to Thomas Tram).

* The 'include' directory contains the declaration files (.h) for all the C files in the 'source', 'main' and 'tools' directories.

* The 'test' directory contains executable programs to test the outputs of SONG.

* The 'scripts' directory contains bash and gnuplot scripts to run SONG iteratively and to plot its results.

* The 'ini' and 'pre' directories contain, respectively, parameter and precision files that can be fed to SONG.

* The 'output' folder, initially empty, will contain the products of SONG computations, usually in the form of text files.



## CONTACT
For any additional support please email sam.lawrence.sl@googlemail.com


## CITATIONS
Writing SONG took more than four years. If you found SONG useful, please acknowledge our effort and that of our collaborators by citing one or more of the following references. If in doubt, please use the first reference.

* G. W. Pettinari, C. Fidler, R. Crittenden, K. Koyama, and D. Wands. "The intrinsic bispectrum of the cosmic microwave background". J. Cosmology Astropart. Phys., 04(2013)003, doi: 10.1088/1475-7516/2013/04/003 [[3]]
```
@article{pettinari:2013a,
       author = {{Pettinari}, G.~W. and {Fidler}, C. and {Crittenden}, R. and {Koyama}, K. and {Wands}, D.},
        title = "{The intrinsic bispectrum of the cosmic microwave background}",
      journal = {J. Cosmology Astropart. Phys.},
         year = 2013,
        month = apr,
       volume = 4,
          eid = {003},
        pages = {3},
          doi = {10.1088/1475-7516/2013/04/003},
       adsurl = {http://adsabs.harvard.edu/abs/2013JCAP...04..003P},
archivePrefix = "arXiv",
       eprint = {1302.0832}
}
```

* G. W. Pettinari, "The intrinsic bispectrum of the cosmic microwave background [PhD Thesis]". Springer Theses Series, Volume XXIII, 2015 [[10]]
```
@book{pettinari:2015a,
	   Author = {Pettinari, G.~W.},
	Publisher = {Springer International Publishing},
	   Series = {Springer Theses},
	    Title = {The Intrinsic Bispectrum of the Cosmic Microwave Background},
	      Url = {http://www.springer.com/gp/book/9783319218816},
	     Year = {2015},
          doi = {10.1088/1475-7516/2013/04/003},
         Isbn = {9783319218823},
archivePrefix = "arXiv",
       eprint = {1405.2280}
}
```

* G. W. Pettinari, C. Fidler, R. Crittenden, K. Koyama, A. Lewis & D. Wands. "Impact of polarization on the intrinsic cosmic microwave background bispectrum". Phys. Rev. D., 90, 103010, doi: 10.1103/PhysRevD.90.103010 [[4]]
```
@article{pettinari:2014a,
	   author = {Pettinari, G.~W. and Fidler, C. and Crittenden, R. and Koyama, K. and Lewis, A and Wands, D.},
	    issue = {10},
	  journal = {Phys. Rev. D},
	    month = {Nov},
	 numpages = {6},
	    pages = {103010},
	publisher = {American Physical Society},
	    title = {Impact of polarization on the intrinsic cosmic microwave background bispectrum},
	      url = {http://link.aps.org/doi/10.1103/PhysRevD.90.103010},
	   volume = {90},
	     year = {2014},
archivePrefix = "arXiv",
       eprint = {1406.2981},
          doi = {10.1103/PhysRevD.90.103010}
}
```


[1]: http://camb.info/ "The CAMB Boltzmann code"
[2]: http://class-code.net/ "The CLASS Boltzmann code"
[3]: http://arxiv.org/abs/1302.0832 "The intrinsic bispectrum of the Cosmic Microwave Background"
[4]: http://arxiv.org/abs/1406.2981 "Impact of polarisation on the intrinsic CMB bispectrum"
[5]: http://arxiv.org/abs/1312.4448 "Spectral distortions in the cosmic microwave background polarization"
[6]: http://arxiv.org/abs/1401.3296 "The intrinsic B-mode polarisation of the Cosmic Microwave Background"
[7]: http://www2.iap.fr/users/pitrou/cmbquick.htm "The CMBquick 2nd-order Boltzmann code"
[8]: http://arxiv.org/abs/1212.3573 "The CMB bispectrum from recombination"
[9]: http://arxiv.org/abs/1104.2933 "The Cosmic Linear Anisotropy Solving System (CLASS)"
[10]: http://arxiv.org/abs/1405.2280 "The intrinsic bispectrum of the Cosmic Microwave Background [PhD thesis]"
[11]: http://arxiv.org/abs/1511.07801 "A precise numerical estimation of the magnetic field generated around recombination"
[Macports]: https://www.macports.org
[Homebrew]: http://brew.sh
[HPC]: http://hpc.sourceforge.net
[Github tutorials]: https://help.github.com/articles/good-resources-for-learning-git-and-github/ "Good Resources for Learning Git and GitHub"
[Github page]: https://github.com/coccoinomane/song.git
[master branch]: https://github.com/coccoinomane/song/tree/master
[develop branch]: https://github.com/coccoinomane/song/tree/develop


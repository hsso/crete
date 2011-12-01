# Make instructions for sky (Solaris/SunOS)

# (c) Michiel Hogerheijde / Floris van der Tak 2000
#     michiel@strw.leidenuniv.nl, vdtak@sron.rug.nl
#     http://www.strw.leidenuniv.nl/~michiel
#     http://www.sron.rug.nl/~vdtak
#
#     This file is part of the 'ratran' molecular excitation and
#     radiative transfer code. The one-dimensional version of this code
#     is publicly available; the two-dimensional version is available on
#     collaborative basis. Although the code has been thoroughly
#     tested, the authors do not claim it is free of errors or that it
#     gives correct results in all situations. Any publication 
#     making use of this code should include a reference to 
#     Hogerheijde & van der Tak, 2000, A&A, 362, 697.

.SILENT:

MIROBJ = $(RATRAN)/sky/mmain.f \
      $(RATRAN)/sky/openim.f \
      $(RATRAN)/sky/closeim.f

FITSOBJ = $(RATRAN)/sky/fmain.f \
      $(RATRAN)/sky/openfits.f 

COMMON = $(RATRAN)/sky/getinputs.f \
      $(RATRAN)/sky/losintegr_$(DIM)d.f \
      $(RATRAN)/sky/readmodel.f \
      $(RATRAN)/sky/vfunc_$(DIM)d.f \
      $(RATRAN)/sky/numerical.f \
      $(RATRAN)/sky/numrep2slatec.f \
      $(RATRAN)/sky/slatec_routines.f \
      $(VELO).f \
      $(KAPPA).f

FITSLIB = -L$(CFITSIO) -lcfitsio -L/usr/X/lib -lX11 -lsocket -lnsl -lm

MIRLIBS = `mirlibs`

OPT = -I. -silent -O       

fits.exe: skycommon.inc 
	f77 $(OPT) $(FITSOBJ) $(COMMON) -o $@ $(FITSLIB)
	strip $@

mir.exe: skycommon.inc 
	f77 $(OPT) $(MIROBJ) $(COMMON) -o $@ $(MIRLIBS)
	strip $@


# Make instructions for amc (Solaris, SunOS)

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
#     gives correct results in all situations. Any publication making 
#     use of this code should include a reference to 
#     Hogerheijde & van der Tak, 2000, A&A, 362, 697.

.SILENT:

# compiler options for Solaris (f77)
OPT = -I. -silent -O       

OBJ = $(RATRAN)/amc/main.f \
      $(RATRAN)/amc/blowpops.f \
      $(RATRAN)/amc/getinputs.f \
      $(RATRAN)/amc/getjbar.f \
      $(RATRAN)/amc/getmatrix.f \
      $(RATRAN)/amc/molinit.f \
      $(RATRAN)/amc/numerical.f \
      $(RATRAN)/amc/numrep2slatec.f \
      $(RATRAN)/amc/slatec_routines.f \
      $(RATRAN)/amc/readmodel.f \
      $(RATRAN)/amc/photon_$(DIM)d.f \
      $(RATRAN)/amc/stateq.f \
      $(RATRAN)/amc/vfunc_$(DIM)d.f \
      $(VELO).f \
      $(KAPPA).f

$(LOCAL).exe: amccommon.inc 
	f77 $(OPT) $(OBJ) -o $@
	strip $@

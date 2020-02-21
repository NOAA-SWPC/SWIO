lib_LIBRARIES = libswio.a

libswio_a_SOURCES = swio.F90 swio_data.F90 swio_methods.F90
libswio_a_FCFLAGS =

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libswio_a-swio.$(OBJEXT) : libswio_a-swio_data.$(OBJEXT) libswio_a-swio_methods.$(OBJEXT)

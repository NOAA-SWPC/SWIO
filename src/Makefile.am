SUBDIRS = .

lib_LIBRARIES = libswio.a

libswio_a_SOURCES = swio.F90

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

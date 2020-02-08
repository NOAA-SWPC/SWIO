SUBDIRS = src
dist_doc_DATA = README

ACLOCAL_AMFLAGS = -I m4

nodist_include_HEADERS = src/swio.$(FC_MODEXT)
nodist_data_DATA = swio.mk
swio.mk:
	@echo "# ESMF self-describing build dependency makefile fragment" > $(@)
	@echo >> $(@)
	@echo "ESMF_DEP_FRONT     = swio"                          >> $(@)
	@echo "ESMF_DEP_INCPATH   = $(includedir)"                 >> $(@)
	@echo "ESMF_DEP_CMPL_OBJS = "                              >> $(@)
	@echo "ESMF_DEP_LINK_OBJS = $(libdir)/libswio.a"           >> $(@)
	@echo "ESMF_DEP_SHRD_PATH = $(SWIO_DEP_SHRD_PATH)"         >> $(@)
	@echo "ESMF_DEP_SHRD_LIBS = $(SWIO_DEP_SHRD_LIBS)"         >> $(@)
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" $(@)
	@echo

uninstall-local:
	-rm -f $(datarootdir)/swio.mk

CLEANFILES = swio.mk

## Makefile.am
## tRNAscan-SE
##
## Settings for generating Makefile.in

# not a GNU package
# build object files in their source location
AUTOMAKE_OPTIONS = foreign subdir-objects

# set flags
AM_CFLAGS = -O
LDADD = -lm

# macros for program dependencies
EUFIND = src/eufind_const.h src/eufind_main.c src/pavesi.c
SQUIDUTILS = src/seqstat_main.c src/reformat_main.c src/shuffle_main.c src/revcomp_main.c src/getopt.c
SQUIDSRC = src/alignio.c src/sqerror.c src/sqio.c src/iupac.c src/msf.c src/revcomp.c src/selex.c src/sre_ctype.c src/sre_math.c src/sre_string.c src/stack.c src/types.c src/gnuregex.c src/interleaved.c src/getopt.c
SQUIDHDRS = src/squid.h src/sqfuncs.h src/gnuregex.h
COVEMAIN = src/align_main.c src/build_main.c src/emit_main.c src/learnprior_main.c src/train_main.c src/scan_main.c src/score_main.c src/structcheck_main.c
COVESRC = src/align.c src/dbviterbi.c src/debug.c src/emit.c src/fast-dbviterbi.c src/fastmodelmaker.c src/konings.c src/lengthdist.c src/maxmodelmaker.c src/misc.c src/model.c src/modelmaking.c src/prior.c src/probify.c src/save.c src/smallviterbi.c src/scorestack.c src/structs.c src/trace.c src/viterbi.c $(SQUIDSRC)

# compile binaries
bin_PROGRAMS = bin/eufindtRNA bin/trnascan-1.4 bin/covels-SE bin/coves-SE
bin_eufindtRNA_SOURCES = $(SQUIDHDRS) $(SQUIDSRC) $(EUFIND)
bin_trnascan_1_4_SOURCES = src/trnascan.c
bin_trnascan_1_4_CFLAGS = -DTSCANDIR=\"$(libdir)\/tRNAscan-SE\/models\"  $(AM_CFLAGS)
bin_covels_SE_SOURCES = $(COVESRC) src/scan_main.c
bin_coves_SE_SOURCES = $(COVESRC) src/score_main.c

# process scripts
# treat tRNAscan-SE.conf as a script since it's just easier
dist_bin_SCRIPTS = bin/sstofa bin/fasta2gsi bin/EukHighConfidenceFilter bin/MitoHighConfidenceFilter tRNAscan-SE tRNAscan-SE.conf
# for a few files, need to substitute perl path, bin directory, and lib directory
# first, define a function for it
do_subst = sed -e 's,\@bindir\@,$(bindir),g' -e 's,\@libdir\@,$(libdir),g' -e 's,[@]PERL[@],$(PERL),g' -e 's,\@INFERNAL\@,$(INFERNAL),g' 

# extra steps for cleaning/distribution
CLEANFILES = $(dist_bin_SCRIPTS) $(srcdir)/lib/*
include_HEADERS = $(srcdir)/src/*.h 
EXTRA_DIST = lib src/sstofa.in src/fasta2gsi.in EukHighConfidenceFilter MitoHighConfidenceFilter tRNAscan-SE.src tRNAscan-SE.conf.src 

# then create build rules

tRNAscan-SE: tRNAscan-SE.src Makefile
	$(do_subst) < $(srcdir)/tRNAscan-SE.src > tRNAscan-SE
	chmod +x tRNAscan-SE

tRNAscan-SE.conf: tRNAscan-SE.conf.src Makefile
	$(do_subst) < $(srcdir)/tRNAscan-SE.conf.src > tRNAscan-SE.conf
	chmod +x tRNAscan-SE.conf

bin/sstofa: src/sstofa.in Makefile
	$(do_subst) < $(srcdir)/src/sstofa.in > bin/sstofa
	chmod +x bin/sstofa

bin/fasta2gsi: src/fasta2gsi.in Makefile
	$(do_subst) < $(srcdir)/src/fasta2gsi.in > bin/fasta2gsi
	chmod +x bin/fasta2gsi

bin/EukHighConfidenceFilter: EukHighConfidenceFilter.in Makefile
	$(do_subst) < $(srcdir)/EukHighConfidenceFilter.in > bin/EukHighConfidenceFilter
	chmod +x bin/EukHighConfidenceFilter

bin/MitoHighConfidenceFilter: MitoHighConfidenceFilter.in Makefile
	$(do_subst) < $(srcdir)/MitoHighConfidenceFilter.in > bin/MitoHighConfidenceFilter
	chmod +x bin/MitoHighConfidenceFilter

# copy libraries. We don't have "libraries" in the libgmp.so.4 kind, but we do have data libraries, so we don't have to deal with lib_LTLIBRARIES.
install-data-hook:
	chmod u+w $(exec_prefix)
	test -d $(libdir) || mkdir $(libdir)
	chmod u+w $(libdir)
	test -d $(libdir)/tRNAscan-SE || mkdir $(libdir)/tRNAscan-SE
	cp -rf $(srcdir)/lib/* $(libdir)/tRNAscan-SE 

# remove lib/ files during uninstallation. A bit of a hack, since lib files are write-protected by default. Files 
uninstall-hook:
	chmod -R u+w $(libdir)/tRNAscan-SE
	rm -rf $(libdir)/tRNAscan-SE


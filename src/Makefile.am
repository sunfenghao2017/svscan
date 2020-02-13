bin_PROGRAMS = svscan svtools bedsolp

svscan_LDADD = $(LDFLAGS) -lhts -lz -lpthread -lxlsxwriter -lbwa

svscan_SOURCES = aligner.cpp breakpoint.cpp annotator.cpp dpbamrecord.cpp junction.cpp \
		 stats.cpp bcfreport.cpp fusionopt.cpp main.cpp msa.cpp options.cpp \
		 region.cpp srbamrecord.cpp svrecord.cpp svscanner.cpp tsvreporter.cpp \
		 cgranges.cpp bed.cpp bamtotb.cpp realnfilter.cpp extanno.cpp

svtools_LDADD = $(LDFLAGS) -lhts -lz -lpthread -lxlsxwriter -lbwa

svtools_SOURCES = svtools.cpp svbam.cpp svdnadb.cpp fusewlist.cpp fusionreporter.cpp \
		  svmerge.cpp svfilter.cpp fusionopt.cpp svrnadb.cpp svdebug.cpp \
		  cgranges.cpp svcreg.cpp bamtotb.cpp realnfilter.cpp extanno.cpp \
		  rsvg.cpp

bedsolp_LDADD = $(LDFLAGS) -lz

bedsolp_SOURCES = cgranges.cpp olpbedl.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}

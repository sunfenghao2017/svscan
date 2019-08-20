bin_PROGRAMS = sver svmerge

sver_LDADD = $(LDFLAGS)

sver_SOURCES = aligner.cpp breakpoint.cpp annotator.cpp dpbamrecord.cpp junction.cpp stats.cpp bcfreport.cpp dptra.cpp \
	       main.cpp msa.cpp options.cpp region.cpp srbamrecord.cpp svrecord.cpp svscanner.cpp tsvreporter.cpp

svmerge_LDADD = $(LDFLAGS)

svmerge_SOURCES = mergevcf.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}

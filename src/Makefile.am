bin_PROGRAMS = sver svmerge getAlnByZF prepAnnoDB

sver_LDADD = $(LDFLAGS)

sver_SOURCES = aligner.cpp breakpoint.cpp annotator.cpp dpbamrecord.cpp junction.cpp stats.cpp bcfreport.cpp dptra.cpp \
	       main.cpp msa.cpp options.cpp region.cpp srbamrecord.cpp svrecord.cpp svscanner.cpp tsvreporter.cpp

svmerge_LDADD = $(LDFLAGS)

svmerge_SOURCES = mergevcf.cpp

getAlnByZF_LDADD = $(LDFLAGS)

getAlnByZF_SOURCES = getAlnByZF.cpp

prepAnnoDB_LDADD = $(LDFLAGS)

prepAnnoDB_SOURCES = prepAnnoDB.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}

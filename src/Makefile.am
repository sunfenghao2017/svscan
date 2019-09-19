bin_PROGRAMS = svscan svmerge fusereport svfilter getAlnByZF prepAnnoDB prepFuseWhiteList

svscan_LDADD = $(LDFLAGS)

svscan_SOURCES = aligner.cpp breakpoint.cpp annotator.cpp dpbamrecord.cpp junction.cpp stats.cpp bcfreport.cpp fusionopt.cpp \
	         main.cpp msa.cpp options.cpp region.cpp srbamrecord.cpp svrecord.cpp svscanner.cpp tsvreporter.cpp

svmerge_LDADD = $(LDFLAGS)

svmerge_SOURCES = mergevcf.cpp

fusereport_LDADD = $(LDFLAGS)

fusereport_SOURCES = fusionreporter.cpp fusionopt.cpp

svfilter_LDADD = $(LDFLAGS)

svfilter_SOURCES = svfilter.cpp

getAlnByZF_LDADD = $(LDFLAGS)

getAlnByZF_SOURCES = getAlnByZF.cpp

prepAnnoDB_LDADD = $(LDFLAGS)

prepAnnoDB_SOURCES = prepAnnoDB.cpp

prepFuseWhiteList_LDADD = $(LDFLAGS)

prepFuseWhiteList_SOURCES = prepFuseWhiteList.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}

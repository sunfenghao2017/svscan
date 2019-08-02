bin_PROGRAMS = sver

sver_LDADD = $(LDFLAGS)

sver_SOURCES = aligner.cpp breakpoint.cpp annotator.cpp dpbamrecord.cpp junction.cpp stats.cpp bcfreport.cpp \
	       main.cpp msa.cpp options.cpp region.cpp srbamrecord.cpp svrecord.cpp svscanner.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}

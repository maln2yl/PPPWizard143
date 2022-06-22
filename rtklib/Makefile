CC=gcc
CFLAGS=-O3 -DENAGLO -DENAGAL -DNFREQ=4 -DNEXOBS=1

#****************

# Objects
FICLIB = rtklib.a
FICOBJ = \
	binex.o\
	convkml.o\
	convrnx.o\
	crescent.o\
	datum.o\
	download.o\
	ephemeris.o\
	geoid.o\
	gw10.o\
	ionex.o\
	javad.o\
	lambda.o\
	novatel.o\
	nvs.o\
	options.o\
	pntpos.o\
	postpos.o\
	ppp.o\
	ppp_ar.o\
	preceph.o\
	qzslex.o\
	rcvlex.o\
	rcvraw.o\
	rinex.o\
	rtcm.o\
	rtcm2.o\
	rtcm3.o\
	rtcm3e.o\
	rtkcmn.o\
	rtkpos.o\
	rtksvr.o\
	sbas.o\
	skytraq.o\
	solution.o\
	ss2.o\
	stream.o\
	streamsvr.o\
	tle.o\
	ublox.o\
	septentrio.o\
	atmo.o

#################
# Compilation

$(FICLIB) : $(FICOBJ)
	@rm -f $@ 
	@echo "Creating rtklib.a..."
	@ar -ru $@ $(FICOBJ)
	@ranlib $@

%.o : %.c
	@rm -f $*.o
	@$(CC) -c $(CFLAGS) $^
    
clean :
	rm -f *.o $(FICLIB)


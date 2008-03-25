COMP=g++

BINDIR = /usr/local/bin

CPROGS = Surface.cpp
LIBS =

PROGS = ${CPROGS:.cpp=} 
OPROGS = ${CPROGS:.cpp=.o}

all: $(PROGS) 

install:
	-cp $(PROGS)  $(BINDIR)

clean:
	rm -f *.o core ${PROGS}

Surface : Surface.cpp Surface.o  
	$(COMP) $(C_FLAGS) -o Surface Face.cpp FastLoop.cpp Loop.cpp PatchDisk.cpp Segment.cpp Surface.cpp Vertex.cpp utils.cpp $(LIBS)

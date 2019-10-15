CC = g++
CFLAGS = -std=c++11 -fPIC
LFLAGS = -L/home/wondrash/Dropbox/c++/lib
INCLUDES = -I/home/wondrash/Dropbox/c++/include
LIBS = -lprops
STATLIB = /home/wondrash/c++/lib/libprops.a
SRCS = main.cpp factory.cpp base_classes.cpp residues.cpp
OBJS = $(SRCS:.cpp=.o)

default: confgen 
	@echo "all binaries have been compiled"

confgen: $(OBJS)
	$(CC) $(INCLUDES) $(CFLAGS) $(OBJS) $(STATLIB) -static -o confgen

%.o:%.cpp
	$(CC) $(INCLUDES) $(CFLAGS) -c $< -o $@
	
clean:
	rm confgen
	rm $(OBJS);

#ViennaRNA package location
ViennaRNA=~/software/include/ViennaRNA
ViennaLIB=~/software/lib/libRNA.a

NAME=graphSim

OBJ = $(NAME)_cmdline.o\
			main.o\
			graph.o

# should stay unchanged: 
CPP = g++
CC = gcc
CPPFLAGS =  -O2 -g -Wall -std=c++0x -fexceptions -Wno-write-strings
CFLAGS =  -O2 -g -Wall -fexceptions -Wno-write-strings
LFLAGS = -O2 -fopenmp

DIRS = -I $(ViennaRNA)

LIBS = $(ViennaLIB)

all: $(OBJ)
	$(CPP) $(LFLAGS) $(DIRS) $(OBJ) $(LIBS) -o $(NAME)
	rm -f $(OBJ)

$(NAME)_cmdline.h $(NAME)_cmdline.c: $(NAME).ggo
	gengetopt -i $(NAME).ggo

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(DIRS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(DIRS) -c $<

clean:
	rm -f $(OBJ)
	rm -f $(NAME)
	rm -f $(NAME)_cmdline.c $(NAME)_cmdline.h


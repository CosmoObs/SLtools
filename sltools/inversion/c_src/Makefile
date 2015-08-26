MAIN = main.cpp

SOURCES = create_grid.cpp \
          pseudo_elliptical_models.cpp \
          radial_models.cpp \
          composit.cpp
OBJ = $(SOURCES:.cpp=.o)


CC = g++
RM = rm

CFLAGS = -Wall -Weffc++ -Wextra -std=c++0x -lpthread

gbclib:
	$(CC) $(CFLAGS) -c $(SOURCES)

clean:
	$(RM) $(OBJ)


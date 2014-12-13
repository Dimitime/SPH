# Makefile for the project
CC = g++
LIBS = -lGL -lGLU -lGLEW -lglfw3 -lm -lXrandr -lXi -lX11 -lXxf86vm -lpthread

all: sph

sph: main.cpp
	$(CC) -o sph main.cpp -Wall $(LIBS)

clean:
	rm -rf *o sph

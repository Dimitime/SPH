# Makefile for the project
CC = g++ -std=c++0x
INCLUDE = -I./deps/
LIBS = -lGL -lGLU -lGLEW -lglut

all: sph

sph: main.cpp
	$(CC)  -o sph main.cpp SphUtils.cpp -Wall $(INCLUDE) $(LIBS) 

clean:
	rm -rf *o sph

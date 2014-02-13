all : CURVE.o nurbs.o
	g++ -o CURVE.o nurbs.o `pkg-config --cflags opencv --libs` -lGLU -lglut -lGL
nurbs.o : ../source/nurbs.cpp ../source/nurbs.h
	g++ -o nurbs.o -c ../source/nurbs.cpp 
SURFACE.o : SURFACE.cpp ../source/nurbs.h
	g++ -o SURFACE.o -c SURFACE.cpp 
CURVE.o : CRUVE_FIT.cpp ../sourcenurbs.h
	g++ -o CURVE.o -c CURVE_FIT.cpp
clean :
	rm *.o
#Way to generate static linking library
libnurbs.a : nurbs.o
	ar crv libmyhello.a hello.o  
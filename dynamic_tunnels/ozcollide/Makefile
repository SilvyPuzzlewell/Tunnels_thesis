
#CXX=ccache g++ -ggdb3 -Wall
CXX=ccache g++ -O2
CXXINCLUDE+=-I/usr/local/include -I../
LDFLAGS+=-L/usr/local/lib -lpthread -lm
LIBTARGET=libozcollide.a

SRC=$(wildcard *.cpp)
OBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

sources = aabbtree.cpp \
          aabbtree_aabb.cpp \
          aabbtree_poly.cpp \
          aabbtree_sphere.cpp \
          aabbtreeaabb_builder.cpp \
          aabbtreepoly_builder.cpp \
          aabbtreesphere_builder.cpp \
          box.cpp \
          dataio.cpp \
          dist_pointline.cpp \
          ellipsoid.cpp \
          frustum.cpp \
          intr_boxbox.cpp \
          intr_frustumsphere.cpp \
          intr_linebox.cpp \
          intr_lineline.cpp \
          intr_segmenttri.cpp \
          intr_spherebox.cpp \
          intr_spheretri.cpp \
          intr_sphereline.cpp \
          intr_tribox.cpp \
          intr_tripoint.cpp \
          matrix.cpp \
          obb.cpp \
          ozcollide.cpp \
          plane.cpp \
          polygon.cpp \
          sphere.cpp \
          vec2f.cpp \
          vec3f.cpp \
          vector.cpp \
          dist_pointbox.cpp 

objss=$(patsubst %.cpp,%.o,$(sources))


all: $(OBJS)
	ar cru $(LIBTARGET) $(OBJS) 
	chmod +x $(LIBTARGET)


.cpp.o:
	$(CXX) $(CXXINCLUDE) -c $<


clean:
	rm -f $(OBJS) $(LIBTARGET);



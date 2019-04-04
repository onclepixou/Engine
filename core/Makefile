CC=g++

SRCDIR  = src
OBJDIR  = obj
BINDIR  = bin

LINCDIR = inc
INCDIR= -I $(LINCDIR) -I /usr/include -I /usr/local/include


CFLAGS=  -std=c++11 -g -fPIC $(INCDIR)

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(LINCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)



$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	$(CC) $(CFLAGS) -c $< -o $@ -shared


tennin.so: $(OBJECTS)
	$(CC) $(CFLAGS)  -o $(BINDIR)/$@ $^ -shared


.PHONY: clean

clean :
	rm $(OBJDIR)/*.o

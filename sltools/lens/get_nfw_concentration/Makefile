TARGET = get_nfw_concentration
INTERFACE = get_nfw_concentration.i
WRAPFILE = $(TARGET)_wrap.c
WRAPOBJ = $(TARGET)_wrap.o
MODULE = _$(TARGET).so
PYMODULE = $(TARGET).py
SRCS = get_nfw_concentration.c rand.c
COBJS = $(SRCS:.c=.o)
SWIGOPT = -python
PYTHONBIN = $(shell echo `which python`)
PYTHONINC = $(shell echo "import sys; print sys.prefix+'/include/python'+str(sys.version_info[0])+'.'+str(sys.version_info[1])" | ${PYTHONBIN})

all: get_nfw_concentration

get_nfw_concentration: rand.c get_nfw_concentration.c
	gcc -Wall -lm $(SRCS) -o $(TARGET)
module_lnx: $(SRCS) $(INTERFACE)
	swig $(SWIGOPT) $(INTERFACE)
	gcc -Wall -fpic -c $(SRCS) $(WRAPFILE) -I$(PYTHONINC) 
	gcc -lm -shared $(COBJS) $(WRAPOBJ) -o $(MODULE)

module_mac: $(SRCS) $(INTERFACE)
	swig $(SWIGOPT) $(INTERFACE)
	gcc -I$(PYTHONINC) -framework Python -bundle -bundle_loader $(PYTHONBIN) $(SRCS) $(WRAPFILE) -o $(MODULE)

clean: 
	@rm -f $(TARGET) $(WRAPFILE) $(WRAPOBJ) $(COBJS) $(MODULE) $(PYMODULE)


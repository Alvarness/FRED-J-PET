#==================================================#
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
# do something Linux-y
PLUGINLIB=-L. -Wl,--whole-archive -lFredXfuncInterface_linux -Wl,--no-whole-archive 
endif
ifeq ($(UNAME), Darwin)
# do something Darwin-y
PLUGINLIB =-L. -Wl,-all_load -lFredXfuncInterface_mac
endif
#==================================================#


plugin-lib:
	g++ -O2 -Wall -shared -fPIC -o libFredXfunc.so myplugin.cpp $(PLUGINLIB)

clean:
	rm -f *.o libFredXfunc.so

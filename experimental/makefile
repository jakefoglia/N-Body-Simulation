Makefile: 

#PROJECT_DIR = /cygdrive/d/Git/N-Body-Simulation/experimental
#OBJ_DIR = $(PROJECT_DIR)/objs

platform = linux

ifeq ($(platform), linux) # Linux parameters
	CC = g++
	CC_FLAGS = -pthread -O3 #-Dgtest_disable_pthreads=ON


#	LIBS = 

#	INCLUDES =
#	LIB_INCLUDES = 

#	INC = $(INCLUDES) $(LIBINCLUDES)  

	BUILD = sim


else					# Windows Parameters
	CC = g++
	CC_FLAGS = -O3 -Dgtest_disable_pthreads=ON

	BUILD = sim.exe

endif





OBJS = main.o

main.o: main.cpp 
	$(CC) $(CC_FLAGS) -c $^
	mv $@ objs/$@

clean: 
	rm -f $(OBJ_DIR)/*.o	

sim: clean $(OBJS)
	$(CC) $(CC_FLAGS) -o $(BUILD) objs/$(OBJS)


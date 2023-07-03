WARN = -Wall -Wsign-compare -Wextra
# add -Warray-bounds=2 for run-time memory checks (does it actually work?)
# OPT  = -O2 -g -pg # for use with gprof
## EXTRA may be set on the command line, e.g., "EXTRA=''" or "EXTRA='-pg -DNDEBUG'"
# EXTRA = -DNDEBUG -Wno-unused-variable -Wno-unused-parameter  
# OPT  = -g -march=native -mtune=native -O3 ${PRODUCTION}
EXTRA=-DNDEBUG
# use -DNDEBUG to disable assert, -Wfatal-errors to stop at 1st error
OPT = -g -march=native -mtune=native -O3 ${EXTRA}
CC = gcc
CFLAGS = ${OPT} ${WARN}
CPP = g++

src = vecdec.c utils_m4ri.c  utils.c ${VECDEC_HDRS}

#-Wmissing-prototypes   -Wconversion -Wextra
# -std=gnu11
VECDEC_HDRS = vecdec.h util_m4ri.h utils.h mmio.h
VECDEC_SECS = vecdec.c util_m4ri.c utils.c mmio.c
VECDEC_OBJS = vecdec.o util_m4ri.o utils.o mmio.o
## unused (documentation only): mainpage.h

src = Makefile ${VECDEC_HDRS} ${VECDEC_SRC} 

default: give_help


# The main decoder program
vecdec: ${VECDEC_OBJS} ${VECDEC_HDRS} Makefile 
	${CPP} ${VECDEC_OBJS} -o vecdec -l m4ri

# depend on all headers 
${VECDEC_OBJS}: ${VECDEC_HDRS} Makefile

mmio.o: mmio.c Makefile
	${CC} ${OPT} -c $<

clean:
	rm -f *~ *.tmp *.out *.o ../input/*~ TAGS

dox: 
	doxygen ../Doxyfile mainpage.h

tags:
	etags *.c *.h 


veryclean: clean
	rm -f vecdec vecdec.exe vecdec.zip 

zip: ${src}
	zip -9 -o vecdec.zip ${src} ../Doxyfile ../README.md

give_help:
	@echo -e "\e[36m make vecdec \e[0m" "   to compile \`vecdec\` \e[32m a simple vectorized decoder\e[0m"
	@echo -e "\e[36m make clean \e[0m   " "to\e[31m remove\e[0m object files and such"
	@echo -e "\e[36m make veryclean \e[0m" "also\e[31m remove\e[0m executables and zip file"
	@echo -e "\e[36m make zip \e[0m    " "to make \e[35m \`vecdec.zip\`\e[0m file with C/C++ sources"
	@echo -e "\e[36m make dox \e[0m    " "to run doxygen to prepare \e[35m documentation\e[0m in 'docs/' \e[31m(broken)\e[0m"
	@echo -e "\e[36m make tags \e[0m   " "to create \e[35m TAGS\e[0m file"


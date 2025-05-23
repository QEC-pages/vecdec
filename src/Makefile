WARN = -Wall -Wsign-compare -Wextra -Wimplicit-fallthrough -Wno-expansion-to-defined
# add -Warray-bounds=2 for run-time memory checks (does it actually work?)
# OPT  = -O2 -g -pg # for use with gprof
# use `-ggdb3` for use with valgrind
## EXTRA may be set on the command line, e.g., "EXTRA=''" or "EXTRA='-pg -DNDEBUG'"
# EXTRA = -DNDEBUG -Wno-unused-variable -Wno-unused-parameter  
# OPT  = -g -march=native -mtune=native -O3 
EXTRA=-DNDEBUG
VER=-DUSE_QLLR
# use -DNDEBUG to disable assert, -Wfatal-errors to stop at 1st error

# define this on the command line for any other compilation options 
# ADD = 

# to compile on Mac you may need to remove '-march=native'
OPT = -g -march=native -mtune=native -O3 ${ADD} ${EXTRA}

# if you have M4RI library outside normal search path
# uncomment the line below (you may edit the location if necessary)
# INC = -I../../m4ri -L../../m4ri/.libs
#
CFLAGS = ${OPT} ${VER} ${WARN} ${INC} 
CC = gcc
CPP = g++


#-Wmissing-prototypes   -Wconversion -Wextra
# -std=gnu11
VECDEC_HDRS= vecdec.h util_m4ri.h utils.h mmio.h tinymt64.h uthash.h qllr.h vec.h
VECDEC_SECS= vecdec.c util_m4ri.c utils.c mmio.c tinymt64.c qllr.c dec_iter.c dec_ris.c dec_pre.c star_poly.c 
VECDEC_OBJS= vecdec.o util_m4ri.o utils.o tinymt64.o qllr.o dec_iter.o dec_ris.o dec_pre.o star_poly.o
VECDEC_XTRA= mmio.o 
## unused (documentation only): mainpage.h

src = Makefile ${VECDEC_HDRS} ${VECDEC_SRC} 

default: give_help

all: vecdec mtx_sub mtx_qc mtx_rep mtx_fw mtx_01_to_nz

# The vecdec program
vecdec: ${VECDEC_OBJS} ${VECDEC_XTRA} ${VECDEC_HDRS} Makefile 
	${CPP} ${CFLAGS} ${VECDEC_OBJS} ${VECDEC_XTRA} -o vecdec -l m4ri

%.o: %.c ${VECDEC_HDRS} Makefile 
	${CC} ${CFLAGS} -c $<

mtx_sub: mtx_sub.c utils.h util_m4ri.h util_m4ri.o utils.o mmio.o Makefile 
	${CC} ${OPT} ${WARN} -o mtx_sub $< util_m4ri.o utils.o mmio.o -l m4ri -lm

mtx_qc: mtx_qc.c utils.h util_m4ri.h util_m4ri.o utils.o mmio.o Makefile 
	${CC} ${OPT} ${WARN} -o mtx_qc $< util_m4ri.o utils.o mmio.o -l m4ri -lm

mtx_rep: mtx_rep.c utils.h util_m4ri.h util_m4ri.o utils.o mmio.o Makefile 
	${CC} ${OPT} ${WARN} -o mtx_rep $< util_m4ri.o utils.o mmio.o -l m4ri -lm

mtx_fw: mtx_fw.c utils.h vec.h Makefile tinymt64.o util_m4ri.o utils.o mmio.o
	${CC} ${OPT} ${WARN} -o mtx_fw $< tinymt64.o util_m4ri.o utils.o mmio.o -l m4ri -l m 

mtx_01_to_nz: mtx_01_to_nz.c Makefile 
	${CC} ${OPT} ${WARN} -o mtx_01_to_nz $< 

clean:
	rm -f *~ *.tmp *.out *.o ../input/*~ TAGS

dox: 
	doxygen ../Doxyfile mainpage.h

tags:
	etags *.c *.h 


veryclean: clean
	rm -f vecdec vecdec.exe vecdec.zip vecdec_dbl

zip: ${src}
	zip -9 -o vecdec.zip ${src} ../Doxyfile ../README.md

give_help:
	@echo  "\e[36m make -j all \e[0m" "   to compile all programs in the collection"
	@echo  "\e[36m make vecdec \e[0m" "   to compile \`vecdec\` \e[32m a simple vectorized decoder\e[0m"
	@echo  "\e[36m make clean \e[0m   " "to\e[31m remove\e[0m object files and such"
	@echo  "\e[36m make veryclean \e[0m" "also\e[31m remove\e[0m executables and zip file"
	@echo  "\e[36m make zip \e[0m    " "to make \e[35m \`vecdec.zip\`\e[0m file with C/C++ sources"
	@echo  "\e[36m make dox \e[0m    " "to run doxygen to prepare \e[35m documentation\e[0m in 'docs/' \e[31m(broken)\e[0m"
	@echo  "\e[36m make tags \e[0m   " "to create \e[35m TAGS\e[0m file"
	@echo
	@echo  "\e[31m The program does not compile?\e[0m"
	@echo  "  - Make sure you have '\e[32mlibm4ri_dev\e[0m' on your system."
	@echo  "    Please see '../README.md' for installation instructions."
	@echo  "  - On a Mac, you may want to disable '\e[32m-march=native\e[0m' compiler switch by running" 
	@echo  "\e[36m make vecdec OPT='-g -mtune=native -O3'\e[0m" 


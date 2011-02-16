#
# Build all routines
#


CODE_DIRS       = `ls | grep -v DLM`

default : all

all clean linux irix osf:
	@ for d in ${CODE_DIRS}; do \
		(if [ -d $$d ]; \
			then cd $$d; \
			echo; \
			echo Making $$d; \
			${MAKE} $@; \
		fi); \
	done

# The old way where we listed all directories 
#all clean linux irix osf:
#	@ for d in ${CODE_DIRS}; do \
#		(cd $$d; echo;echo Making $$d; ${MAKE} $@); \
#	done

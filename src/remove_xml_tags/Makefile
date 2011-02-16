include make.conf

CC	= gcc
LD	= gcc
INCLUDE	= -I${IDL_DIR}/external/include
CCFLAGS = ${INCLUDE} ${M_CFLAGS}
LDFLAGS = ${M_LDFLAGS}

LIBOBJ  = remove_xml_tags.o
LIBRARY = remove_xml_tags.so
DLM     = remove_xml_tags.dlm

all: ${LIBRARY}
	cp ${LIBRARY} ../DLM
	cp ${DLM} ../DLM
clean:
	rm -f ${LIBRARY} ${LIBOBJ}
	rm -f ../DLM/${LIBRARY}
	rm ../DLM/${DLM}

${LIBRARY}: ${LIBOBJ}
	${LD} ${LDFLAGS} -o $@ $?

.c.o:
	${CC} -c ${CCFLAGS} $<

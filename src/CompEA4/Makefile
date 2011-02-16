#  Makefile for total_int.  This file is part of sdssidl
#  Copyright (C) 2005  Erin Sheldon, NYU.  erin.sheldon at gmail.com
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# This creates the M_CFLAGS and M_LDFLAGS
include make.conf

CC	= gcc
LD	= gcc
INCLUDE	= -I${IDL_DIR}/external/include
CCFLAGS = ${M_CFLAGS} ${INCLUDE}
LDFLAGS = ${M_LDFLAGS}

LIBOBJS  = CompEA4.o CompEA4IDL.o
LIBRARY = CompEA4IDL.so
DLM      = CompEA4IDL.dlm

all: ${LIBRARY}
	cp ${LIBRARY} ../DLM/
	cp ${DLM} ../DLM/
clean:
	- rm -f ${LIBRARY} ${LIBOBJS}
	- rm -f ../DLM/${LIBRARY}
	- rm -f ../DLM/${DLM}

${LIBRARY}: ${LIBOBJS}
	@ echo linking $@
	${LD} ${LDFLAGS} -o $@ ${LIBOBJS}

.c.o:
	@ echo compiling $<
	${CC} -c ${CCFLAGS} $<

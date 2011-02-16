Below is a file containing the clusters. There are two cavets:
	I made no attempt to search for AGN, so the x-ray flux may or may
		not be coming from the cluster
	I made a cut on position: the brightest cluster galaxy had to be
		within 4 sigma of the RASS position, using the RASS position
		error. Clearly this is too conservative. A better approach
		might be to ask if the RASS source is within a 1/2 Mpc of the
		brightest cluster galaxy, but that is not what I did for
		this sample.

The file has:
         name:		rass concated with ra,dec of the RASS source
         ra:		ra of the brightest cluster galaxy
         dec:		dec of the brightest cluster galaxy
         z:		photo-z of cluster
         run:		run of BCG candidate
         camCol:	camera column of BCG candidate
         field:		field of BCG candidate
         row:		image row of BCG candidate
         col:		image col of BCG candidate
         bcgRa:		id number of BCG candidate
         bcgDec:	junk	(duplicate of dec)


===================

typedef struct {
	 char name[40];
	 char ra[40];
	 char dec[40];
	 float z;
	 int run;
	 float camCol;
	 float field;
	 float row;
	 float col;
	 float bcgRa;
	 float bcgDec;
} CLUSTER;

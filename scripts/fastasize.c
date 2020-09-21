/*
 * Copyright (c) 2003-2006 University of Chicago and Fellowship
 * for Interpretations of Genomes. All Rights Reserved.
 *
 * This file is part of the SEED Toolkit.
 * 
 * The SEED Toolkit is free software. You can redistribute
 * it and/or modify it under the terms of the SEED Toolkit
 * Public License. 
 *
 * You should have received a copy of the SEED Toolkit Public License
 * along with this program; if not write to the University of Chicago
 * at info@ci.uchicago.edu or the Fellowship for Interpretation of
 * Genomes at veronika@thefig.info or download a copy from
 * http://www.theseed.org/LICENSE.TXT.
 */


/*  fastasize.c
 *
 *  Usage:  fastasize -t < fasta_file  > nseq nresidues
 *  or      fastasize    < fasta_file  > id nresidues ...
 */
/*  These include files are appropriate for Machintosh OS X  */

#include <stdio.h>
#include <ctype.h>   /*  isspace() */
#include <stdlib.h>  /*  exit()    */
#include <unistd.h>  /*  read()    */

#define  BUFLEN    (256*1024)
#define  INPLEN    ( 64*1024)
#define  IDLEN     ( 16*1024)
#define  DFLT_INDEX_INTERVAL  10000

#define  fillbuf(buf, len)  read( (int) 0, (void *) buf, (size_t) len )

/*  Function prototypes:  */

void report_seq( char * id, unsigned long seqlen );

void report_ttl( int n_seq, unsigned long long ttllen );

void usage( char *prog );


unsigned char  buffer[BUFLEN];
char           idbuf[IDLEN+1];


int main ( int argc, char **argv ) {
    unsigned long long  ttllen;
    unsigned char  *bptr;
    unsigned long   c, seqlen;
    int             totalonly, n_seq, idlen, ntogo;

    /* -t flag returns only the total */

    totalonly = 0;
    if ( ( argc > 1 ) && ( argv[1][0] == '-'  )
                      && ( argv[1][1] == 't'  )
                      && ( argv[1][2] == '\0' )
       ) {
	totalonly = 1;
	argc--;
    }
    else if (argc > 1 ) { usage( argv[0] ); }

    idbuf[0] = '\0';  /* initialize to empty string */
    bptr   = buffer;  /* pointer to next character in buffer */
    ntogo  =  0;      /* unused characters in input buffer */
    n_seq  =  0;
    seqlen =  0;
    ttllen =  0;

    /*  Process one line of input  */

    while ( 1 ) {
	if ( ntogo <= 0 ) {
	    if ( ( ntogo = fillbuf( buffer, BUFLEN ) ) <= 0 ) {
		if ( totalonly ) {
		    ttllen += seqlen;
		    report_ttl( n_seq, ttllen );
		}
		else {
		    report_seq( idbuf, seqlen );
		}
		exit( ntogo );
	    }
	    bptr= buffer;
	}
	c = *bptr++; ntogo--;

	/*  Line could start with >, or be sequence data  */

	if ( c == '>' ) {

	    /*  New sequence.  Is there an previous sequence to report?  */

	    if ( ! totalonly ) report_seq( idbuf, seqlen );

	    /*  Adjust cumulative values and reset sequence values  */

	    ttllen += seqlen;
	    seqlen = 0;

	    if ( ntogo <= 0 ) {
		if ( ( ntogo = fillbuf( buffer, BUFLEN ) ) <= 0 ) {
		    if ( totalonly ) report_ttl( n_seq, ttllen );
		    exit( ntogo );
		}
		bptr = buffer;
	    }
	    c = *bptr++; ntogo--;

	    /*  Make a copy of the new id  */

	    idlen = 0;
	    while ( ( ! isspace(c) ) && ( idlen < IDLEN ) ) {
		idbuf[ idlen++ ] = c;
		if ( ntogo <= 0 ) {
		    if ( ( ntogo = fillbuf( buffer, BUFLEN ) ) <= 0 ) {
		        if ( totalonly ) report_ttl( n_seq, ttllen );
		        exit(ntogo);
		    }
		    bptr = buffer;
		}
		c = *bptr++; ntogo--;
	    }
	    idbuf[ idlen ] = '\0';

	    /*  report truncated id  */

	    if ( ! isspace(c) ) {
		fprintf( stderr, "Sequence id truncated to %d characters:\n", (int) IDLEN );
		fprintf( stderr, ">%s\n", idbuf );
	    }

	    /*  Flush the rest of the input line  */

	    while ( c != '\n' ) {
		if ( ntogo <= 0 ) {
		    if ( ( ntogo = fillbuf( buffer, BUFLEN ) ) <= 0 ) {
		        if ( totalonly ) report_ttl( n_seq, ttllen );
		        exit(ntogo);
		    }
		    bptr = buffer;
		}
		c = *bptr++; ntogo--;
	    }

	    n_seq++;  /*  First data for a new sequence  */
	}

	/*  Not an id line, so it's data:  */

	else {
	    while ( c != '\n' ) {  /* finish the line */
		if ( ! isspace( c ) ) {
		    seqlen++;
		}

		/*  Next character  */

		if ( ntogo <= 0 ) {
		    if ( ( ntogo = fillbuf( buffer, BUFLEN ) ) <= 0 ) {
		        if ( totalonly ) {
		            ttllen += seqlen;
		            report_ttl( n_seq, ttllen );
		        }
			else {
			    report_seq( idbuf, seqlen );
			}
		        exit( ntogo );
		    }
		    bptr = buffer;
		}
		c = *bptr++; ntogo--;
	    }
	}
	/*  Go back to top to process next line  */
    }

    exit( 0 );   /*  never get here  */
}


void report_seq( char * id, unsigned long seqlen ) {
    if ( id && id[0] ) {
	printf( "%s\t%lu\n", id, seqlen );
    }
}


void report_ttl( int n_seq, unsigned long long ttllen ) {
    printf( "%d\t%llu\n", n_seq, ttllen );
}


void usage( char *prog ) {
    fprintf( stderr,
            "Usage:  %s -t < fasta_file  > nseq \\t nresidues\n"
            "or      %s    < fasta_file  > id \\t nresidues\\n ...\n",
             prog, prog
           );
    exit(0);
}

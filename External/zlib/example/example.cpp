#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "zlib/zlib/zlib.h"

#define LENGTH 100000000 // memory block size for reading

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
		exit(-1);
	}

	unsigned char *buffer = (unsigned char *) malloc(LENGTH * sizeof(*buffer));

	if (! buffer) {
		fprintf (stderr, "Allocation of %d bytes failed: %s.\n", LENGTH, strerror(errno));
		exit (EXIT_FAILURE);
    }

	gzFile fh = gzopen(argv[1], "rb");
	if (! fh) {
        fprintf (stderr, "gzopen of '%s' for reading failed: %s.\n", argv[1], strerror(errno));
        exit (EXIT_FAILURE);
    }

	gzFile oh = gzopen(argv[2], "wbT");
	if (! oh) {
        fprintf (stderr, "gzopen of '%s' for uncompressed writing failed: %s.\n", argv[2], strerror(errno));
        exit (EXIT_FAILURE);
    }
	
	unsigned int x = 0xf7142a0e;
    while (1) {
        int err;                    
        int bytes_read;
        bytes_read = gzread(fh, buffer, LENGTH);
        gzwrite(oh, buffer, sizeof(*buffer) * bytes_read);
		for (int i = 0; i < bytes_read; i++) x = (x << 8) ^ buffer[i];
		fprintf(stderr, "Read %d bytes.\n", bytes_read); fflush(stderr);
        if (bytes_read < LENGTH) {
            if (gzeof (fh)) {
				fprintf(stderr, "EOF reached\n");
                break;
            }
            else {
				fprintf(stderr, "Incompatible partial read and EOF status\n");
                const char * error_string;
                error_string = gzerror (fh, & err);
                if (err) {
                    fprintf (stderr, "Error: %s.\n", error_string);
                    exit (EXIT_FAILURE);
                }
            }
        }
    }
    gzclose(fh);
	gzclose(oh);
	fprintf(stderr, "Signature: %0x\n", x);

    return 0;
}


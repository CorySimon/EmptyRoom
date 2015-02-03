/*
 * write_to_outputfile.h
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */

#ifndef WRITE_TO_OUTPUTFILE_H_
#define WRITE_TO_OUTPUTFILE_H_

void write_settings_to_outputfile(FILE * outputfile, GridParameters parameters, Framework framework) {
    //
	// Write time
	//
	time_t t = time(0); // get time now
    struct tm * now = localtime ( & t );
    fprintf(outputfile, "date: %d-%d-%d", (now->tm_year + 1900), (now->tm_mon + 1), now->tm_mday);
}


#endif /* WRITE_TO_OUTPUTFILE_H_ */

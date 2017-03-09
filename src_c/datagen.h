/*
 * datagen.h
 *
 *  Created on: Feb 22, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DATAGEN_H_
#define SRC_C_DATAGEN_H_

#include "source.h"

// This is the main routine. The other sub functions are used to suppor DataGen.
void DataGen();

// Inspiral source functions
void Print_Source(source_t* source);
void Load_Source(source_t* source);


#endif /* SRC_C_DATAGEN_H_ */

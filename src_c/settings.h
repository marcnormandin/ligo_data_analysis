/*
 * settings.h
 *
 *  Created on: Mar 29, 2017
 *      Author: marcnormandin
 */

#ifndef TESTS_SETTINGS_FILE_SETTINGS_H_
#define TESTS_SETTINGS_FILE_SETTINGS_H_

#include <stdio.h>

#define SETTING_MAX_KEY_SIZE 255
#define SETTING_MAX_VAL_SIZE 255

typedef struct setting_s {
	char key[SETTING_MAX_KEY_SIZE];
	char val[SETTING_MAX_VAL_SIZE];

	/* linked list */
	struct setting_s *next;

} setting_t;

typedef struct settings_file_s {
	FILE *fid;
	char *fname;

	setting_t *first;

} settings_file_t;


settings_file_t* settings_file_open(char *filename);
void settings_file_close(settings_file_t* sf);
char* settings_file_get_value(settings_file_t *sf, char *key);
void settings_file_print(settings_file_t *sf);

#endif /* TESTS_SETTINGS_FILE_SETTINGS_H_ */

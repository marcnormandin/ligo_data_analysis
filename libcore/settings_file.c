#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "settings_file.h"

setting_t* read_setting(settings_file_t* sf) {
	assert(sf != NULL);

	char key[SETTING_MAX_KEY_SIZE];
	char val[SETTING_MAX_VAL_SIZE];
	memset(key, '\0', SETTING_MAX_KEY_SIZE *sizeof(char));
	memset(val, '\0', SETTING_MAX_VAL_SIZE *sizeof(char));

	setting_t *s;
	int r;

	r = fscanf(sf->fid, "%s %s", key, val);

	if (r == EOF) {
		return NULL;
	}

	if (r != 2) {
		printf("Error: Invalid formatted setting line. Must be of the form [key] [val]\n");
		return NULL;
	}

	s = (setting_t*) malloc(sizeof(setting_t));
	if (s == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for setting in read_setting(). Exiting.\n");
		exit(-1);
	}

	strncpy(s->key, key, SETTING_MAX_KEY_SIZE);
	strncpy(s->val, val, SETTING_MAX_VAL_SIZE);
	s->next = NULL;

	return s;
}

void add_setting(settings_file_t *sf, setting_t *s) {
	assert(sf != NULL);
	assert(s != NULL);

	setting_t *current;
	current = sf->first;
	if (current == NULL) {
		sf->first = s;
		return;
	}

	while (current->next != NULL) {
		current = current->next;
	}

	current->next = s;
}

void read_file(settings_file_t* sf) {
	assert(sf != NULL);

	setting_t *s;
	assert(sf);

	while ((s = read_setting(sf)) != NULL) {
		add_setting(sf, s);
	}
}

settings_file_t* settings_file_open(const char *filename) {
	assert(filename != NULL);

	FILE *fid;
	settings_file_t *sf;

	assert(filename);

	fid = fopen(filename, "r");
	if (fid == NULL) {
		fprintf(stderr, "Error: Unable to open the settings file (%s) for reading. Exiting.\n", filename);
		exit(-1);
	}

	sf = (settings_file_t*) malloc(sizeof(settings_file_t));
	if (sf == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for settings_file_t. Exiting.\n");
		exit(-1);
	}
	sf->fid = fid;
	sf->first = NULL;

	/* fixme */
	size_t n = strnlen(filename, 254);
	sf->fname = (char*) malloc( (n+1) * sizeof(char));
	memset(sf->fname, '\0', (n+1) * sizeof(char));
	strncpy(sf->fname, filename, n);

	read_file(sf);

	return sf;
}

void settings_file_close(settings_file_t* sf) {
	assert(sf != NULL);

	setting_t *current, *next;

	if (sf->fid != NULL) {
		current = sf->first;
		while(current != NULL) {
			next = current->next;
			free(current);
			current = next;
		}
		fclose(sf->fid);
	}

	assert(sf->fname != NULL);
	free(sf->fname);
	sf->fname = NULL;

	free(sf);
}

const char* settings_file_get_value(settings_file_t *sf, const char *key) {
	assert(sf != NULL);
	assert(key != NULL);

	setting_t *current;
	current = sf->first;
	while (current != NULL) {
		if (strcmp(key, current->key) == 0) {
			return current->val;
		}
		current = current->next;
	}
	return NULL;
}

void settings_file_print(settings_file_t *sf) {
	assert(sf != NULL);

	setting_t *current;
	current = sf->first;
	while (current != NULL) {
		printf("%s: %s\n", current->key, current->val);
		current = current->next;
	}
}

int settings_file_num_settings(settings_file_t *sf) {
	assert(sf != NULL);

	setting_t *current;
	int count;

	assert(sf != NULL);

	count = 0;

	current = sf->first;

	while(current != NULL) {
		count++;
		current = current->next;
	}

	return count;
}

const char* settings_file_get_key_by_index(settings_file_t *sf, size_t index) {
	assert(sf != NULL);

	setting_t *current;
	int i;

	assert(sf);

	i = 0;
	current = sf->first;
	while (current != NULL) {
		if (i == index) {
			return current->key;
		}
		current = current->next;
		i++;
	}

	fprintf(stderr, "Error: Attempt to return an invalid index. Exiting.\n");
	exit(-1);
}

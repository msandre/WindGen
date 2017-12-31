/* Wind simulation library.
 * Copyright (C) 2013  Michael Andre.
 * Copyright (C) 2013  Technische Universitaet Muenchen.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "input.h"
#include <regex.h>

/**
 * read_file:
 * Read a file into a buffer and return a pointer to the first element.
 */
char * read_file(const char * file, char comment)
{
  FILE * fp;
  int fpos, pos, c, size;
  char * buf;
  fp = fopen(file,"r");
  if (fp == NULL) { return NULL; }

  size = 0;
  while(fgetc(fp) != EOF)
    size = size + 1;

  rewind(fp);

  buf = (char *) malloc(size+1);
  pos = 0;
  for (fpos=1; fpos <= size; fpos++) {
    c = fgetc(fp);

    if (c == (int) comment) { // jump to next line                                                                         
      while (c != '\n' && fpos < size) {
        c = fgetc(fp);
        fpos = fpos + 1;
      }

      if (fpos < size) {
        c = fgetc(fp);
        fpos = fpos + 1;
      }
    }  

    buf[pos] = (char) c;
    pos = pos + 1;
  }

  buf[pos] = '\0';
  fclose(fp);
  return buf;
}


/**
 * reg_match:
 * Match a regular expression in a string. If a match is found,
 * copy it into a buffer and return a pointer to the first element.
 * Otherwise, return NULL.
 */
char * reg_match(const char * str, const char * pattern)
{
  int stat, pos, size;
  char * buf;
  regex_t re;
  regmatch_t match[1];

  stat = regcomp(&re, pattern, REG_EXTENDED);
  if (stat != 0) { return NULL; }

  stat = regexec(&re, str, 1, match, 0);
  if (stat != 0) { return NULL; }

  size = match[0].rm_eo - match[0].rm_so;
  buf = (char *) malloc( size+1 );

  for (pos=0; pos < size; pos++) {
    buf[pos] = str[match[0].rm_so + pos];
  }

  buf[pos] = '\0';
  regfree(&re);
  return buf;
}


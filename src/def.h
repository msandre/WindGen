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

/*
 * Author   : michael.andre@tum.de
 * Created  : 2013-06-12
 */

#ifndef __DEF_H__
#define __DEF_H__

#define WFS_ERROR_SIZE 256
#define WFS_ERROR(err, msg) { snprintf(err, WFS_ERROR_SIZE, msg); return NULL; }

#endif /* __DEF_H__ */

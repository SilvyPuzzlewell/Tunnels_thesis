/*
OZCollide - Collision Detection Library
Copyright (C) 2006-2010 by Igor Kravtchenko

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

Contact the author: igor@tsarevitch.org
*/

#ifndef OZCOLLIDE_H
#define OZCOLLIDE_H

#ifdef OZCOLLIDE_DLLEXPORT
 #define OZCOLLIDE_API __declspec(dllexport)
 #pragma message("OZCollide - DLL Export")
#elif OZCOLLIDE_DLLIMPORT
 #define OZCOLLIDE_API __declspec(dllimport)
 #pragma message("OZCollide - DLL Import")
#else
 #define OZCOLLIDE_API
#endif

#ifndef ENTER_NAMESPACE_OZCOLLIDE
#define ENTER_NAMESPACE_OZCOLLIDE namespace ozcollide {
#endif

#ifndef LEAVE_NAMESPACE
#define LEAVE_NAMESPACE }
#endif

#define ozinline __inline

#define OZCOLLISON_VERSION 0x0100

#ifndef SAFE_ACOS
#define SAFE_ACOS(x) ( (x) >= 1 ? 0 : ( (x) <= -1 ? OZ_PI : acos(x) ) )
#endif
#ifndef SAFE_ASIN
#define SAFE_ASIN(x) ( (x) >= 1 ? OZ_HALFPI : ( (x) <= -1 ? -OZ_HALFPI : asin(x) ) )
#endif

#define RAD_TO_DEG(x) ((x) * 180 / OZ_PI)
#define DEG_TO_RAD(x) ((x) * OZ_PI / 180)

#define OZ_PI			3.14159265359f
#define OZ_HALFPI		(OZ_PI * 0.5)
#define OZ_QUARTPI		(OZ_PI * 0.25)
#define OZ_TWOPI		(OZ_PI * 2)

#define OZ_COS45	0.707106781187f
#define OZ_SQR2		1.41421356237f
#define OZ_SQR3		1.73205080757f
#define OZ_INVQSQR2	0.707106781187f
#define OZ_INVQSQR3	0.57735026919f

#define BI_MID(a,b,c,d)		( ((a)<<24) + ((b)<<16) + ((c)<<8) + (d) )
#define LI_MID(a,b,c,d)		( ((d)<<24) + ((c)<<16) + ((b)<<8) + (a) )
#define MID					LI_MID

#define SAFE_FREE(p) if (p) { ::free(p); p = NULL; }
#define SAFE_DELETE(p) if (p) { delete p; p = NULL; }
#define SAFE_ARRAYDELETE(p) if (p) { delete [] p; p = NULL; }

ENTER_NAMESPACE_OZCOLLIDE

/**
 * Determines an error occured during an operation.
 */
enum ERR {
    NOERR = 0, /**< No error.*/
    ERR_CANNOT_OPEN = 0x11, /**< Cannot open the requested ressource.*/
    ERR_INVALID_FORMAT = 0x12, /**< The provided data are not in the correct format.*/
};

OZCOLLIDE_API const char* getErrorString(ERR);
OZCOLLIDE_API int getLibVersion();

LEAVE_NAMESPACE

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aabbtree.h"
#include "aabbtree_aabb.h"
#include "aabbtree_poly.h"
#include "aabbtree_sphere.h"
#include "aabbtreeaabb_builder.h"
#include "aabbtreepoly_builder.h"
#include "aabbtreesphere_builder.h"
#include "box.h"
#include "dataio.h"
#include "dist_pointbox.h"
#include "dist_pointline.h"
#include "ellipsoid.h"
#include "frustum.h"
#include "intr_boxbox.h"
#include "intr_frustumsphere.h"
#include "intr_linebox.h"
#include "intr_lineline.h"
#include "intr_segmenttri.h"
#include "intr_spherebox.h"
#include "intr_sphereline.h"
#include "intr_spheretri.h"
#include "intr_tribox.h"
#include "intr_tripoint.h"
#include "matrix.h"
#include "monitor.h"
#include "obb.h"
#include "plane.h"
#include "polygon.h"
#include "sphere.h"
#include "vec2f.h"
#include "vec3f.h"
#include "vector.h"

#endif

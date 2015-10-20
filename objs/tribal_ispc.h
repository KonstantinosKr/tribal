//
// objs/tribal_ispc.h
// (Header automatically generated by the ispc compiler.)
// DO NOT EDIT THIS FILE.
//

#ifndef ISPC_OBJS_TRIBAL_ISPC_H
#define ISPC_OBJS_TRIBAL_ISPC_H

#include <stdint.h>



#ifdef __cplusplus
namespace ispc { /* namespace */
#endif // __cplusplus

///////////////////////////////////////////////////////////////////////////
// Functions exported from ispc code
///////////////////////////////////////////////////////////////////////////
#if defined(__cplusplus) && !defined(__ISPC_NO_EXTERN_C)
extern "C" {
#endif // __cplusplus
    extern void generate_triangles_and_velocities(double * lo, double * hi, int64_t nt, double * t[][3], double *  * v, uint64_t * tid, uint64_t * pid);
    extern void generate_velocities(double * lo, double * hi, int64_t nt, double *  * v);
    extern void integrate_triangles(double step, double * lo, double * hi, int64_t nt, double * t[][3], double *  * v);
    extern void scale_triangles(double * lo, double * hi, int64_t nt, double * t[][3]);
#if defined(__cplusplus) && !defined(__ISPC_NO_EXTERN_C)
} /* end extern C */
#endif // __cplusplus


#ifdef __cplusplus
} /* namespace */
#endif // __cplusplus

#endif // ISPC_OBJS_TRIBAL_ISPC_H

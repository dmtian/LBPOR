/*
 * API for SHAKE implementation.
 */

#ifndef SHAKE_H__
#define SHAKE_H__

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================== */
/*
 * SHAKE implementation.
 *
 * The computation uses a caller-provided context structure. Call sequence:
 *
 *   - Initialise the structure with shake_init(). The "capacity" defines
 *     the variant which is used; it must be a multiple of 64, at least 64,
 *     at most 1600, but the extreme values don't make much sense.
 *     SHAKE128 uses capacity 256; SHAKE256 uses capacity 512.
 *
 *   - Inject some data with shake_inject(). This can be done with an
 *     arbitrary number of calls; each call injects an arbitrary number
 *     of bytes.
 *
 *   - Call shake_flip() to switch to extraction mode.
 *
 *   - Obtain successive output bytes with shake_extract(). Each call
 *     may ask for an arbitrary number of bytes (even zero).
 *
 * There is no explicit deallocation, since there is no dynamic
 * allocation.
 *
 * Caller is responsible for ensuring that a context is properly
 * initialised, that shake_flip() is called, that no call to
 * shake_inject() is made after shake_flip (and before the next
 * shake_init() on the same structure), and that no call to
 * shake_extract() is made before shake_flip().
 */

/*
 * SHAKE context. Contents are opaque (don't access them directly).
 */
typedef struct {
	unsigned char dbuf[200];
	size_t dptr;
	size_t rate;
	uint64_t A[25];
} shake_context;

/*
 * Initialise (or reinitialise) the provided context structure, for a
 * given capacity. Capacity is expressed in bits, and is typically
 * twice the "bit security level" you are aiming at. In practice,
 * use 256 for SHAKE128, 512 for SHAKE256.
 */
void shake_init(shake_context *sc, int capacity);

/*
 * Inject some data bytes in the provided SHAKE context (which must have
 * been initialised previously with shake_init(), and not yet "flipped").
 */
void shake_inject(shake_context *sc, const void *data, size_t len);

/*
 * "Flip" the provided SHAKE context to extraction mode. The context must
 * have been previously initialised, and not yet flipped.
 */
void shake_flip(shake_context *sc);

/*
 * Extract bytes from a flipped SHAKE context. Several successive calls
 * can be made, yielding the next bytes of the conceptually infinite
 * output.
 */
void shake_extract(shake_context *sc, void *out, size_t len);

/* ==================================================================== */

#ifdef __cplusplus
}
#endif

#endif

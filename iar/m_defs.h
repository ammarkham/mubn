/* 
 * File:   m_defs.h
 * Author: Andres
 *
 * Created on April 12, 2011, 12:26 AM
 */

#ifndef M_DEFS_H
#define	M_DEFS_H

#ifdef	__cplusplus
extern "C" {
#endif

//MSP430 has a 16-bit word size
#define WORDSIZE 16
// If m = ceil(log_2 (p)) is the bit length of p and t = ceil(m/W) is the wordlength.
// For now I'll restrict to work with numbers in [0,p-1] where p = p_192 with wordlength 12
#define WORDLENGTH 12

//#define	DEBUG_PRINTF

#ifdef	__cplusplus
}
#endif

#endif	/* M_DEFS_H */


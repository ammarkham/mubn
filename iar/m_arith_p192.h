/* 
 * File:   m_arith_p192.h
 * Author: amolina
 *
 * Created on April 12, 2011, 3:53 PM
 */

#ifndef M_ARITH_P192_H
#define	M_ARITH_P192_H

#ifdef	__cplusplus
extern "C" {
#endif


/**
 * Fast reduction modulo p = p_192 = 2^{192} - 2^{64} -1 c.f Alg. 2.27
 * Input: c, s.t 0<= c <= p^2
 * Output: c mod p
 */
void reduce_mod_p(uint16_t * c, uint16_t * p, uint16_t * c_out);

void multiply_mod_p_192(uint16_t * c, uint16_t * a, uint16_t * b);

/**
 * Implementation of the left to right modular exponentiation algorithm
 * as described in the HAC book by Menezes et.al.
 *
 * @param A The result of raising g to the power of e
 * @param g an element of Z*_p
 * @param e a single precission exponent
 */
void mod_exp_p_192_lr(uint16_t * A, uint16_t * g, uint16_t e);

/**
 * Implementation of the left to right modular exponentiation algorithm
 * as described in the HAC book by Menezes et.al.
 *
 * @param A The result of raising g to the power of e
 * @param g an element of Z*_p
 * @param e a multi-precission exponent
 * @param e_legth the wordlength of the multi-precission exponent
 */
void mod_exp_p_192(uint16_t * A, uint16_t * g, uint16_t * e, uint16_t e_length);

#ifdef	__cplusplus
}
#endif

#endif	/* M_ARITH_P192_H */


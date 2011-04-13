/* 
 * File:   m_arith.h
 * Author: Andres
 *
 * Created on April 12, 2011, 12:27 AM
 */

#ifndef M_ARITH_H
#define	M_ARITH_H

#ifdef	__cplusplus
extern "C" {
#endif
#include <stdint.h>


// A big number will be stored as an array A of WORDLENGTH  WORDSIZE-words where A[0] is the
//least significant bit.

/**
 * One word addition with a carry bit
 * c_i = (a_1 + b_i + epsilon_prime) mod 0xFFFF
 */
uint16_t add_word(uint16_t * c_i, uint16_t a_i, uint16_t b_i, uint16_t epsilon_prime);
/**
 * One word subtraction with a borrow bit
 * c_i = (a_1 - b_i - epsilon_prime) mod 0xFFFF
 */
uint8_t subtract_word(uint16_t * c_i, uint16_t a_i, uint16_t b_i, uint8_t epsilon_prime);
/**
 * Multiprecision addition c.f. Alg. 2.5
 * Input: a, b \in [0,2^{Wt})
 * Output: (epsilon, c) where c = a + b mod 2^{Wt} and epsilon is the carry bit
 */
uint16_t add_mp_elements(uint16_t * pfe_c, uint16_t * pfe_a, uint16_t * pfe_b, uint8_t wordlength);
/**
 * Multiprecision subtraction c.f. Alg. 2.6
 * Input: a, b \in [0,2^{Wt})
 * Output: (epsilon, c) where c = a - b mod 2^{Wt} and epsilon is the borrow bit
 */
uint16_t subtract_mp_elements(uint16_t * pfe_c, uint16_t * pfe_a, uint16_t * pfe_b, uint8_t wordlength);
/**
 * Addition in F_p c.f. Alg. 2.7
 * Input: a, b \in [0,p-1)
 * Output: c = a + b mod p
 */
void add_mod_p(uint16_t * c, uint16_t * a, uint16_t * b, uint16_t * p, uint8_t wordlength);
/**
 * subtraction in F_p c.f. Alg. 2.8
 * Input: a, b \in [0,p-1)
 * Output: c = a - b mod p
 */
void subtract_mod_p(uint16_t * c, uint16_t * a, uint16_t * b, uint16_t * p, uint8_t wordlength);
/**
 * Multiply two single words into a double word
 * Input: a,b words
 * Output: uv one 2-word
 */
void multiply_words(uint16_t a_i, uint16_t b_i, uint16_t * uv);
/**
 * Multiprecision Multiplication c.f. Alg. 2.9
 * Input: a, b \in [0,p-1]
 * Output: c = a*b
 */
void multiply_mp_elements(uint16_t * c, uint16_t * a, uint16_t * b, uint8_t wordlength);
/**
 * Multiprecision Multiplication c.f. Alg. 2.9
 * Input: a, b \in [0,p-1]
 * Output: c = a*b
 */
void multiply_mp_elements2(uint16_t * c, uint16_t * a, uint8_t wordlength_a, uint16_t * b, uint8_t wordlength_b) ;
/**
 * Sets a bn to zero
 */
void set_to_zero(uint16_t * c, uint8_t wordlength);


/**
 * Compares two big nums
 * Returns 1 if a >= b 0 otherwise
 */
uint8_t compare_mp_elements(uint16_t * a, uint16_t * b, uint8_t wordlength);

int are_mp_equal(uint16_t * a, uint16_t * b, uint8_t wordlength);

void mult_by_power_of_b(uint16_t * out, uint16_t wordlength_out, uint16_t * a, uint16_t wordlength_a, uint16_t k);

/*
Dividing mod a power of the radix is done by simply shifting the string right
the corresponding number of places

def div_mod_16_power(z,k):
    if (len(z)/4) <= k:
        z_div = "0"
    else:
        z_div = z[0:-4*k]
    return z_div
*/
void div_by_power_of_b(uint16_t * out, uint16_t * a, uint16_t k, uint16_t wordlength);


void multiply_sp_by_mp_element(uint16_t * c, uint16_t a, uint16_t * b, uint16_t wordlength_b);

void copy_mp(uint16_t * out, uint16_t * in, int wordlength);

int ith_bit(uint16_t e, int i);

void mod_pow_of_b(uint16_t * out, uint16_t wordlength_out, uint16_t * a, uint16_t wordlength_a, uint16_t k);

int bit_length(uint16_t e);

int mp_bit_length(uint16_t * e, uint16_t wordlength);

void multiply_mod_p(uint16_t * c, uint16_t * a, uint16_t * b, uint16_t * p, uint8_t wordlength);

int mp_ith_bit(uint16_t * e, int i);

#ifdef	__cplusplus
}
#endif

#endif	/* M_ARITH_H */

